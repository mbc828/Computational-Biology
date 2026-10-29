# Public source snapshot: replace path placeholders before running.
#!/bin/bash
# GSE235929 总控流水线：下载(4并行续传,无超时) -> 全部齐后逐个 salmon quant -> DONE 标记
# 用 nohup 后台运行，日志写 pipeline.log；可安全重复启动（内部步骤均幂等/断点续传）
set -u
cd /mnt/m/shxdpji/gse235929
LOCK=/tmp/gse235929_pipeline.lock
if [ -f "$LOCK" ] && pgrep -F "$LOCK" >/dev/null 2>&1; then
  echo "已在运行，退出"; exit 0
fi
echo $$ > "$LOCK"

echo "=== [PIPELINE START] $(date) ==="

# ---- 阶段 1：下载直到 12/12 ----
while true; do
  missing=0
  while IFS=$'\t' read -r gsm srr group obj url bytes; do
    [ "$gsm" = "gsm" ] && continue
    f="fastq/${srr}.fastq.gz"
    sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
    [ "$sz" = "$bytes" ] || missing=$((missing+1))
  done < fastq_manifest.tsv
  echo "[下载] 缺 $missing/12 $(date)"
  [ "$missing" -eq 0 ] && break
  n=0; pids=()
  while IFS=$'\t' read -r gsm srr group obj url bytes; do
    [ "$gsm" = "gsm" ] && continue
    f="fastq/${srr}.fastq.gz"
    sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
    [ "$sz" = "$bytes" ] && continue
    curl -sL -C - -o "$f" "$url" &
    pids+=($!); n=$((n+1))
    [ $n -ge 4 ] && break
  done < fastq_manifest.tsv
  [ ${#pids[@]} -gt 0 ] && wait "${pids[@]}"
done
echo "=== [下载完成] $(date) ==="

# ---- 阶段 2：gzip 完整性校验 ----
bad=0
while IFS=$'\t' read -r gsm srr group obj url bytes; do
  [ "$gsm" = "gsm" ] && continue
  f="fastq/${srr}.fastq.gz"
  if ! gzip -t "$f" 2>/dev/null; then
    echo "[校验失败] $f —— 重新下载"
    rm -f "$f"; bad=1
  fi
done < fastq_manifest.tsv
if [ $bad -eq 1 ]; then
  echo "存在损坏文件，重启流水线补下"
  exec bash "$0"
fi
echo "=== [gzip 校验通过] $(date) ==="

# ---- 阶段 3：salmon 定量 ----
while IFS=$'\t' read -r gsm srr group obj url bytes; do
  [ "$gsm" = "gsm" ] && continue
  f="fastq/${srr}.fastq.gz"
  out="quant/${srr}"
  [ -f "${out}/quant.sf" ] && { echo "skip ${srr}"; continue; }
  echo "== quant ${srr} (${gsm}, ${group}) $(date) =="
  salmon quant -i ref/salmon_index -l A -r "${f}" \
      --validateMappings --threads 10 \
      -o "${out}" >> pipeline_quant.log 2>&1 \
    && echo "OK ${srr} $(date)" || echo "FAIL ${srr}"
done < fastq_manifest.tsv

echo "=== [PIPELINE DONE] $(date) ===" > PIPELINE_DONE.txt
rm -f "$LOCK"
