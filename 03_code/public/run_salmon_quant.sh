# Public source snapshot: replace path placeholders before running.
#!/bin/bash
# salmon 联合定量：12 个单端 FASTQ × 人+金葡菌双参考索引；可断点重入
set -u
cd /mnt/m/shxdpji/gse235929
IDX=ref/salmon_index
while IFS=$'\t' read -r gsm srr group obj url bytes; do
  [ "$gsm" = "gsm" ] && continue
  f="fastq/${srr}.fastq.gz"
  out="quant/${srr}"
  if [ -f "${out}/quant.sf" ]; then
    echo "skip ${srr}（已定量）"; continue
  fi
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$sz" != "$bytes" ]; then
    echo "wait ${srr}（FASTQ 未齐 $sz/$bytes）"; continue
  fi
  echo "== quant ${srr} (${gsm}, ${group}) $(date) =="
  salmon quant -i "${IDX}" -l A -r "${f}" \
      --validateMappings --threads 8 \
      -o "${out}" > "${out}.log" 2>&1 \
    && echo "OK ${srr}" || echo "FAIL ${srr}（见 ${out}.log）"
done < fastq_manifest.tsv
echo "[QUANT ROUND DONE] $(date)"
