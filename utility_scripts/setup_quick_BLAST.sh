rm -rf quick_blast
mkdir quick_blast
cp ./PPYR_OGS1.0*.fa.gz quick_blast
cd quick_blast
gunzip ./*.fa.gz
sequenceserver -d . -m 
