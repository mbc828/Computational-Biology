# -*- coding: utf-8 -*-
"""admetlab2_batch.py — 16 成分 × ADMETlab 2.0 全端点批量复测（第二平台交叉验证）
流程：GET 表单页取 CSRF → POST smiles → 解析结果页（数值端点+分类端点+色点+类别）
已下载的 HTML 跳过，可断点续跑。"""
import csv, re, time, sys
from pathlib import Path
import urllib.request, urllib.parse, http.cookiejar

BASE = 'https://admetmesh.scbdd.com'
OUT = Path('admet_xval/results_al2')
OUT.mkdir(parents=True, exist_ok=True)

comps = [(r['Compound'], r['SMILES']) for r in csv.DictReader(open('admet_xval/compounds_16.csv', encoding='utf-8-sig'))]
todo = [c for c in comps if not (OUT / f'{c[0]}.html').exists()]
limit = int(sys.argv[1]) if len(sys.argv) > 1 else 16
print(f'共 {len(comps)} 个，待跑 {len(todo)} 个，本轮上限 {limit}')

def fetch(name, smi):
    cj = http.cookiejar.CookieJar()
    op = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
    op.addheaders = [('User-Agent', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)')]
    html = op.open(BASE + '/service/evaluation/index', timeout=30).read().decode('utf-8', 'ignore')
    tok = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', html).group(1)
    data = urllib.parse.urlencode({'csrfmiddlewaretoken': tok, 'smiles': smi, 'method': '1'}).encode()
    req = urllib.request.Request(BASE + '/service/evaluation/cal', data=data,
            headers={'Referer': BASE + '/service/evaluation/index'})
    body = op.open(req, timeout=180).read().decode('utf-8', 'ignore')
    assert len(body) > 100000, f'结果页异常短: {len(body)}'
    (OUT / f'{name}.html').write_text(body, encoding='utf-8')
    return len(body)

ok = fail = 0
for name, smi in todo[:limit]:
    try:
        n = fetch(name, smi)
        print(f'  OK {name}: {n} B'); ok += 1
    except Exception as e:
        print(f'  FAIL {name}: {type(e).__name__} {e}'); fail += 1
    time.sleep(2.5)
print(f'完成 ok={ok} fail={fail}')
