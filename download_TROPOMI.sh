#!/bin/sh

download_tropomi_data() {
    year=$1
    month=$2
    day=$3

    # Format month and day to ensure two digits
    month=$(printf "%02d" $month)
    day=$(printf "%02d" $day)

    filetest="s3://DIAS/Sentinel-5P/TROPOMI/L2__CO____/$year/$month/$day/"
    echo $filetest

    s3cmd --recursive -c ~/.s3cfg ls $filetest > test.txt
    sed '/NRTI/d' ./test.txt > test0.txt
    sed '/.cdl/d' ./test0.txt > test01.txt
    sed 's#^.\{31\}#s3cmd --recursive --skip-existing -c ~/.s3cfg get #g' test01.txt > test2.txt
    local_path="/nobackupp19/bbyrne1/Test_TROPOMI_download/PRODUCT/$year/$month/$day/"
    sed "s#\$# $local_path#g" test2.txt > test3.sh
    sed -i '1s*^*#!/bin/sh\n*' test3.sh
    chmod u+x test3.sh
    ./test3.sh

    # Delete any .cdl files that were downloaded
    find $local_path -name "*.cdl" -type f -delete
}

# Call the function with the provided arguments
download_tropomi_data "$1" "$2" "$3"
