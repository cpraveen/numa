rm -f *.pdf
pdflatex main.tex

npages=3

for ((i = 1 ; i <= npages ; i++ ))
do 
   echo "------------------------$i------------------------"
   pdfjam main.pdf $i -o p$i.pdf
   pdfbb p$i.pdf
   #pdftocairo -svg p$i.pdf p$i.svg
   pdf2svg p$i.pdf p$i.svg
done
