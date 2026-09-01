#!/bin/bash

cd /mnt/g/B_rapa_work/B_rapa_Heat_MetaAnalysis

# 创建结果文件头
echo "Sample,Total_Reads,Pass_Reads,Q20_Rate,Q30_Rate,GC_Content,Read_Length,Source" > results/QC/qc_metrics_all_60.csv

# 定义60个样本
samples=(
SRR12632445 SRR12632457 SRR12632456 SRR12632437 SRR12632436 SRR12632435
SRR12632407 SRR12632406 SRR12632405 SRR12632430 SRR12632429 SRR12632428
SRR12632446 SRR12632444 SRR12632447 SRR12632404 SRR12632455 SRR12632454
SRR12632453 SRR12632452 SRR12632451 SRR12632424 SRR12632422 SRR12632421
SRR12632443 SRR12632442 SRR12632441 SRR26447755 SRR26447754 SRR26447745
SRR26447751 SRR26447750 SRR26447749 SRR26447741 SRR26447740 SRR26447739
SRR26447744 SRR26447743 SRR26447742 SRR26447748 SRR26447747 SRR26447746
SRR26447738 SRR26447753 SRR26447752 SRR21227071 SRR21227073 SRR21227074
SRR21227072 SRR21227069 SRR21227070 SRR21227075 SRR21227076 SRR21227065
SRR12214241 SRR12214243 SRR12214249 SRR12214245 SRR12214247 SRR12214251
)

echo "开始提取60个样本的质控数据..."

for s in "${samples[@]}"; do
    html="data/counts/${s}_fastp.html"
    json="data/counts/${s}_fastp.json"
    log="data/counts/${s}_fastp.log"
    
    # 判断文件类型
    if [ -f "$json" ]; then
        # ===== JSON 提取 =====
        Rscript -e "
            library(jsonlite)
            d <- fromJSON('$json')
            before <- d\$summary\$before_filtering
            after <- d\$summary\$after_filtering
            cat('$s,')
            cat(before\$total_reads, ',')
            cat(after\$total_reads, ',')
            cat(after\$q20_rate * 100, ',')
            cat(after\$q30_rate * 100, ',')
            cat(after\$gc_content * 100, ',')
            cat(ifelse(is.null(after\$read1_mean_length), 0, after\$read1_mean_length), ',')
            cat('json\n')
        " 2>/dev/null >> results/QC/qc_metrics_all_60.csv
        echo "✅ $s (json)"
        
    elif [ -f "$html" ]; then
        # ===== HTML 提取 =====
        # 提取 After filtering 部分的数据
        total=$(grep -A 10 "After filtering" "$html" | grep "total reads" | head -1 | sed -E 's/.*>([0-9.]+) M.*/\1/')
        if [ -n "$total" ]; then
            total=$(echo "$total * 1000000" | bc)
        else
            total=0
        fi
        
        pass=$(grep -A 10 "Filtering result" "$html" | grep "reads passed filters" | head -1 | sed -E 's/.*>([0-9.]+) M.*/\1/')
        if [ -n "$pass" ]; then
            pass=$(echo "$pass * 1000000" | bc)
        else
            pass=0
        fi
        
        q20=$(grep -A 10 "After filtering" "$html" | grep "Q20" | head -1 | sed -E 's/.*\(([0-9.]+)%\).*/\1/')
        q30=$(grep -A 10 "After filtering" "$html" | grep "Q30" | head -1 | sed -E 's/.*\(([0-9.]+)%\).*/\1/')
        gc=$(grep -A 10 "After filtering" "$html" | grep "GC content" | head -1 | sed -E 's/.*>([0-9.]+)%<.*/\1/')
        
        # 如果上述提取失败，尝试用更宽松的方式
        if [ -z "$q20" ]; then
            q20=$(grep "Q20 bases" "$html" | head -1 | sed -E 's/.*\(([0-9.]+)%\).*/\1/')
        fi
        if [ -z "$q30" ]; then
            q30=$(grep "Q30 bases" "$html" | head -1 | sed -E 's/.*\(([0-9.]+)%\).*/\1/')
        fi
        if [ -z "$gc" ]; then
            gc=$(grep "GC content" "$html" | head -1 | sed -E 's/.*>([0-9.]+)%<.*/\1/')
        fi
        
        echo "$s,$total,$pass,$q20,$q30,$gc,0,html" >> results/QC/qc_metrics_all_60.csv
        echo "✅ $s (html)"
        
    elif [ -f "$log" ]; then
        # ===== LOG 提取 =====
        total=$(grep "total reads:" "$log" | head -1 | sed -E 's/.*: ([0-9,]+).*/\1/' | tr -d ',')
        pass=$(grep "passed filters reads:" "$log" | head -1 | sed -E 's/.*: ([0-9,]+).*/\1/' | tr -d ',')
        q20=$(grep "Q20 rate:" "$log" | head -1 | sed -E 's/.*: ([0-9.]+).*/\1/')
        q30=$(grep "Q30 rate:" "$log" | head -1 | sed -E 's/.*: ([0-9.]+).*/\1/')
        gc=$(grep "GC content:" "$log" | head -1 | sed -E 's/.*: ([0-9.]+).*/\1/')
        
        echo "$s,$total,$pass,$q20,$q30,$gc,0,log" >> results/QC/qc_metrics_all_60.csv
        echo "✅ $s (log)"
        
    else
        echo "$s,NA,NA,NA,NA,NA,NA,missing" >> results/QC/qc_metrics_all_60.csv
        echo "❌ $s (缺失)"
    fi
done

echo ""
echo "========================================="
echo "✅ 完成！共处理60个样本"
echo "结果文件: results/QC/qc_metrics_all_60.csv"
echo "========================================="

# 显示统计
echo ""
echo "按来源统计:"
cut -d',' -f8 results/QC/qc_metrics_all_60.csv | tail -n +2 | sort | uniq -c

