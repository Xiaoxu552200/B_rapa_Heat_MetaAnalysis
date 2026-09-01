#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 fastp 的 HTML/JSON/LOG 文件中提取 Q20、Q30 等质控指标
"""

import os
import json
import re
import csv
from bs4 import BeautifulSoup

# 60个样本
samples = [
    "SRR12632445", "SRR12632457", "SRR12632456", "SRR12632437", "SRR12632436",
    "SRR12632435", "SRR12632407", "SRR12632406", "SRR12632405", "SRR12632430",
    "SRR12632429", "SRR12632428", "SRR12632446", "SRR12632444", "SRR12632447",
    "SRR12632404", "SRR12632455", "SRR12632454", "SRR12632453", "SRR12632452",
    "SRR12632451", "SRR12632424", "SRR12632422", "SRR12632421", "SRR12632443",
    "SRR12632442", "SRR12632441", "SRR26447755", "SRR26447754", "SRR26447745",
    "SRR26447751", "SRR26447750", "SRR26447749", "SRR26447741", "SRR26447740",
    "SRR26447739", "SRR26447744", "SRR26447743", "SRR26447742", "SRR26447748",
    "SRR26447747", "SRR26447746", "SRR26447738", "SRR26447753", "SRR26447752",
    "SRR21227071", "SRR21227073", "SRR21227074", "SRR21227072", "SRR21227069",
    "SRR21227070", "SRR21227075", "SRR21227076", "SRR21227065",
    "SRR12214241", "SRR12214243", "SRR12214249", "SRR12214245", "SRR12214247",
    "SRR12214251"
]

counts_dir = "data/counts"
results = []

def parse_json_file(filepath):
    """解析 JSON 文件"""
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        before = data.get('summary', {}).get('before_filtering', {})
        after = data.get('summary', {}).get('after_filtering', {})
        return {
            'total_reads': before.get('total_reads', 0),
            'pass_reads': after.get('total_reads', 0),
            'q20_rate': after.get('q20_rate', 0) * 100,
            'q30_rate': after.get('q30_rate', 0) * 100,
            'gc_content': after.get('gc_content', 0) * 100,
            'read_length': after.get('read1_mean_length', 0),
            'source': 'json'
        }
    except:
        return None

def parse_html_file(filepath):
    """解析 HTML 文件，提取嵌入的 JSON 数据"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # 方法1: 查找包含 JSON 数据的 script 标签
        soup = BeautifulSoup(content, 'html.parser')
        for script in soup.find_all('script'):
            if script.string:
                text = script.string
                # 查找 JSON 对象
                match = re.search(r'\{[^{}]*"passed_filter"[^{}]*\}', text)
                if match:
                    try:
                        data = json.loads(match.group())
                        before = data.get('summary', {}).get('before_filtering', {})
                        after = data.get('summary', {}).get('after_filtering', {})
                        return {
                            'total_reads': before.get('total_reads', 0),
                            'pass_reads': after.get('total_reads', 0),
                            'q20_rate': after.get('q20_rate', 0) * 100,
                            'q30_rate': after.get('q30_rate', 0) * 100,
                            'gc_content': after.get('gc_content', 0) * 100,
                            'read_length': after.get('read1_mean_length', 0),
                            'source': 'html'
                        }
                    except:
                        pass
        
        # 方法2: 直接从 HTML 中提取表格数据
        # fastp HTML 中通常有表格显示 Q20/Q30
        q20_match = re.search(r'Q20[^0-9]*([0-9.]+)%', content)
        q30_match = re.search(r'Q30[^0-9]*([0-9.]+)%', content)
        gc_match = re.search(r'GC content[^0-9]*([0-9.]+)%', content)
        
        if q20_match and q30_match:
            return {
                'total_reads': 0,
                'pass_reads': 0,
                'q20_rate': float(q20_match.group(1)),
                'q30_rate': float(q30_match.group(1)),
                'gc_content': float(gc_match.group(1)) if gc_match else 0,
                'read_length': 0,
                'source': 'html_table'
            }
        return None
    except Exception as e:
        return None

def parse_log_file(filepath):
    """解析 LOG 文件"""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        total_reads = re.search(r'total reads:\s*([0-9,]+)', content)
        pass_reads = re.search(r'passed filters reads:\s*([0-9,]+)', content)
        q20 = re.search(r'Q20 rate:\s*([0-9.]+)', content)
        q30 = re.search(r'Q30 rate:\s*([0-9.]+)', content)
        gc = re.search(r'GC content:\s*([0-9.]+)', content)
        read_len = re.search(r'length of reads:\s*([0-9.]+)', content)
        
        return {
            'total_reads': int(total_reads.group(1).replace(',', '')) if total_reads else 0,
            'pass_reads': int(pass_reads.group(1).replace(',', '')) if pass_reads else 0,
            'q20_rate': float(q20.group(1)) if q20 else 0,
            'q30_rate': float(q30.group(1)) if q30 else 0,
            'gc_content': float(gc.group(1)) if gc else 0,
            'read_length': float(read_len.group(1)) if read_len else 0,
            'source': 'log'
        }
    except:
        return None

print("========================================")
print("提取 fastp 质控指标 (Python 版本)")
print("========================================\n")

for s in samples:
    html_file = os.path.join(counts_dir, f"{s}_fastp.html")
    json_file = os.path.join(counts_dir, f"{s}_fastp.json")
    log_file = os.path.join(counts_dir, f"{s}_fastp.log")
    
    result = None
    source = "unknown"
    
    if os.path.exists(json_file):
        result = parse_json_file(json_file)
        source = "json" if result else "unknown"
    elif os.path.exists(html_file):
        result = parse_html_file(html_file)
        source = "html" if result else "unknown"
    elif os.path.exists(log_file):
        result = parse_log_file(log_file)
        source = "log" if result else "unknown"
    
    if result and result.get('q20_rate', 0) > 0:
        results.append({
            'Sample': s,
            'Total_Reads': result.get('total_reads', 0),
            'Pass_Reads': result.get('pass_reads', 0),
            'Q20_Rate': round(result.get('q20_rate', 0), 2),
            'Q30_Rate': round(result.get('q30_rate', 0), 2),
            'GC_Content': round(result.get('gc_content', 0), 2),
            'Read_Length': round(result.get('read_length', 0), 2),
            'Source': source
        })
        print(f"  ✅ {s} ({source})")
    else:
        print(f"  ❌ {s}")

# 保存 CSV
output_dir = "results/QC"
os.makedirs(output_dir, exist_ok=True)

with open(os.path.join(output_dir, "qc_metrics_60samples_py.csv"), 'w', newline='') as f:
    if results:
        writer = csv.DictWriter(f, fieldnames=['Sample', 'Total_Reads', 'Pass_Reads', 
                                                'Q20_Rate', 'Q30_Rate', 'GC_Content', 
                                                'Read_Length', 'Source'])
        writer.writeheader()
        writer.writerows(results)

print(f"\n✅ 完成！成功提取 {len(results)} / 60 个样本")
print(f"输出文件: results/QC/qc_metrics_60samples_py.csv")
