far --source test.csfastq --target result_right_nw.csfastq --format csfastq --min-overlap 4 --adapters adapters.csfasta --min-readlength 10 --cut-off 1 --trim-end right > /dev/null

a=`diff correct_result_right_nw.csfastq result_right_nw.csfastq`
echo $a
#b=`diff correct_result_right_nw.csfastq.omitted test.csfastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode csfastq, right, nw"

else
echo "Test 1 OK"
fi

far --source test.csfastq --target result_left_nw.csfastq --format csfastq --min-overlap 4 --adapters adapters.csfasta --min-readlength 10 --cut-off 1 --trim-end left  > /dev/null

a=`diff correct_result_left_nw.csfastq result_left_nw.csfastq`
echo $a
#b=`diff correct_result_left_nw.csfastq.omitted test.csfastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode cdfastq, left, nw"
exit -1
else
echo "Test 2 OK"
fi


far --source test.csfastq --target result_any_nw.csfastq --format csfastq --min-overlap 4 --adapters adapters.csfasta --min-readlength 10 --cut-off 1 --trim-end any  > /dev/null

a=`diff correct_result_any_nw.csfastq result_any_nw.csfastq`
echo $a
#b=`diff correct_result_any_nw.csfastq.omitted test.csfastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode any, left, nw"
exit -1
else
echo "Test 3 OK"
fi

far --source test.csfastq --target result_left_tail_nw.csfastq --format csfastq --min-overlap 4 --adapters adapters.csfasta --min-readlength 10 --cut-off 1 --trim-end left_tail  > /dev/null

a=`diff correct_result_left_tail_nw.csfastq result_left_tail_nw.csfastq`
echo $a
#b=`diff correct_result_left_tail_nw.csfastq.omitted test.csfastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode fasta, left_tail, nw"
exit -1
else
echo "Test 4 OK"
fi


far --source test.csfastq --target result_right_tail_nw.csfastq --format csfastq --min-overlap 4 --adapters adapters.csfasta --min-readlength 10 --cut-off 1 --trim-end right_tail  > /dev/null

a=`diff correct_result_right_tail_nw.csfastq result_right_tail_nw.csfastq`
echo $a
#b=`diff correct_result_right_tail_nw.csfastq.omitted test.csfastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode fasta, right_tail, nw"
exit -1
else
echo "Test 5 OK"
fi
