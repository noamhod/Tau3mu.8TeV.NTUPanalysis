dat=`date +"%d%m%Y"`
fname=backup_${dat}.tgz
tar -czf $fname *.C *.h *.sh SummarizeSyst/*.tex SummarizeSyst/Makefile
echo "backed up version $dat to $fname"
