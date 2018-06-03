rm -rf quick_blast
mkdir quick_blast
cp ./*.fa.gz quick_blast
cd quick_blast
gunzip ./*.fa.gz
sequenceserver -d . -m 
