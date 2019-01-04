genmatrix() {
  echo $(date)":    Column Grep: "$1 >> logs.txt
  echo $(date)":    Full Data Table: "$2 >> logs.txt
  echo $(date)":    BinSize: "$3 >> logs.txt
  echo $(date)":    Outfile: "$4 >> logs.txt
  
  
  header $2 > header.txt
  cut -f 2,4,5 $2 > snpPos.tmp
  sed -i '1 s/^.*$/\tchrs\tpos/' snpPos.tmp
  cut -f $(grep -n "$1" header.txt | cut -d : -f 1 | tr '\n' ',' | sed 's/,$//') $2 > cnVal.tmp
  
  
  end=$3
  start=1
  cnt=1
  while [ $start -le $(header cnVal.tmp  | wc -l) ]
  do
    # Check if the range exceeds the number of columns
    if [ $end -ge $(header cnVal.tmp  | wc -l) ]; then
      end=$(header cnVal.tmp  | wc -l)
    fi
    
    # Subset the range of columns
    echo ${start}" - "${end}
    cut -f${start}-${end} cnVal.tmp > cnVal_${cnt}.tmp
    paste -d '\t' snpPos.tmp cnVal_${cnt}.tmp > ${4}_${cnt}.txt
    rm cnVal_${cnt}.tmp
    
    # Update
    cnt=$(($cnt + 1))
    start=$(($end + 1))
    end=$(( $end + $3 ))
  done
  
  rm snpPos.tmp cnVal.tmp header.txt
}

#genmatrix2 'Log R Ratio' full_data_table.txt 50 GNE_LogR
genmatrix 'B Allele Freq' full_data_table.txt 50 GNE_BAF
genmatrix 'GType' full_data_table.txt 50 GNE_Genotype

