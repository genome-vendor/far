far --source test.fasta --target result_right_nw.fasta --format fasta --min-overlap 4 --adapters adapters.fasta --min-readlength 10 --cut-off 1 --trim-end right > /dev/null

a=`diff correct_result_right_nw.fasta result_right_nw.fasta`
echo $a
b=`diff correct_result_right_nw.omitted test.fasta.omitted`

l1=`expr length "$a"`
l2=`expr length "$b"`

if [ $l1 != 0 ] || [ $l2 != 0 ]; then
echo "error testing mode fasta, right, nw"
exit -1
else
echo "Test 1 OK"
fi

far --source test.fasta --target result_left_nw.fasta --format fasta --min-overlap 4 --adapters adapters.fasta --min-readlength 10 --cut-off 1 --trim-end left  > /dev/null

a=`diff correct_result_left_nw.fasta result_left_nw.fasta`
echo $a
b=`diff correct_result_left_nw.omitted test.fasta.omitted`

l1=`expr length "$a"`
l2=`expr length "$b"`

if [ $l1 != 0 ] || [ $l2 != 0 ]; then
echo "error testing mode fasta, left, nw"
exit -1
else
echo "Test 2 OK"
fi


far --source test.fasta --target result_any_nw.fasta --format fasta --min-overlap 4 --adapters adapters.fasta --min-readlength 10 --cut-off 1 --trim-end any  > /dev/null

a=`diff correct_result_any_nw.fasta result_any_nw.fasta`
echo $a
b=`diff correct_result_any_nw.omitted test.fasta.omitted`

l1=`expr length "$a"`
l2=`expr length "$b"`

if [ $l1 != 0 ] || [ $l2 != 0 ]; then
echo "error testing mode any, left, nw"
exit -1
else
echo "Test 3 OK"
fi

far --source test.fasta --target result_left_tail_nw.fasta --format fasta --min-overlap 4 --adapters adapters.fasta --min-readlength 10 --cut-off 1 --trim-end left_tail  > /dev/null

a=`diff correct_result_left_tail_nw.fasta result_left_tail_nw.fasta`
echo $a
b=`diff correct_result_left_tail_nw.omitted test.fasta.omitted`

l1=`expr length "$a"`
l2=`expr length "$b"`

if [ $l1 != 0 ] || [ $l2 != 0 ]; then
echo "error testing mode fasta, left_tail, nw"
exit -1
else
echo "Test 4 OK"
fi


far --source test.fasta --target result_right_tail_nw.fasta --format fasta --min-overlap 4 --adapters adapters.fasta --min-readlength 10 --cut-off 1 --trim-end right_tail  > /dev/null

a=`diff correct_result_right_tail_nw.fasta result_right_tail_nw.fasta`
echo $a
b=`diff correct_result_right_tail_nw.omitted test.fasta.omitted`

l1=`expr length "$a"`
l2=`expr length "$b"`

if [ $l1 != 0 ] || [ $l2 != 0 ]; then
echo "error testing mode fasta, right_tail, nw"
exit -1
else
echo "Test 5 OK"
fi
