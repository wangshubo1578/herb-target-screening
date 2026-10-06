"""
common.py — 中药单体反向筛选 skill 的公共工具模块
提供：HTTP 会话（重试/超时）、JSON 中间文件 IO、日志、路径管理。
"""
import json
import os
import sys
import time
import logging
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ---------------- 常量 ----------------
TCMSP_TOKEN = "fa6966e547446646375e7c2a176163ab"  # 实测恒定，首页不含 token
HERB_API = "http://herb.ac.cn/chedi/api/"
TCMSP_SEARCH = "https://old.tcmsp-e.com/tcmspsearch.php"
TCMSP_TARGET = "https://old.tcmsp-e.com/target.php"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CLINICALTRIALS_API = "https://clinicaltrials.gov/api/v2/studies"
ALPHAFOLD_FILES = "https://alphafold.ebi.ac.uk/files"

HEADERS = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                  "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"
}

# ---------------- 日志 ----------------
def setup_logger(name="herb-screening"):
    logger = logging.getLogger(name)
    if not logger.handlers:
        h = logging.StreamHandler(sys.stderr)
        h.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(h)
        logger.setLevel(logging.INFO)
    return logger

LOG = setup_logger()

def log(msg, level="info"):
    getattr(LOG, level, LOG.info)(msg)

# ---------------- HTTP ----------------
def get_session(retries=3, backoff=1.5):
    """带重试的 requests.Session"""
    s = requests.Session()
    retry = Retry(total=retries, backoff_factor=backoff,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET", "POST"])
    adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=10)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    s.headers.update(HEADERS)
    return s

SESSION = get_session()

def http_get(url, params=None, timeout=30, **kw):
    r = SESSION.get(url, params=params, timeout=timeout, **kw)
    r.raise_for_status()
    return r

def http_post_json(url, payload, timeout=30, **kw):
    r = SESSION.post(url, json=payload, timeout=timeout,
                     headers={"Content-Type": "application/json"}, **kw)
    r.raise_for_status()
    return r

# ---------------- 路径 / JSON IO ----------------
def get_outdir(path="./tcm_run"):
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p

def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")
    log(f"写出 {path.name}")

def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def step_path(outdir, n, name):
    """生成第 N 步的 JSON 路径，如 02_raw_compounds.json"""
    return Path(outdir) / f"{n:02d}_{name}.json"

# ---------------- 容错 ----------------
def safe_step(fn, *args, default=None, **kw):
    """执行一步；失败时记录并返回 default，不中断整链"""
    try:
        return fn(*args, **kw)
    except Exception as e:
        log(f"步骤 {getattr(fn, '__name__', fn)} 失败: {e}", level="warning")
        return default

def normalize_name(name):
    """化合物名归一并去掉连接符，用于黑名单/去重匹配"""
    import re
    return re.sub(r"[-_ ]+", "", (name or "").lower().strip())
