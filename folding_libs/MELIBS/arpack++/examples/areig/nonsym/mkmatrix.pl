#!/usr/bin/perl

#########################################
#
#  A short little perl script to make a
#  generic matrix of size n
#
#  Usage: mkmatrix.pl n > file.matrix
#
#########################################

# size of the matrix
$n = $ARGV[0];
$sparse = $ARGV[$#ARGV];
if ($sparse == 0)
{
  $fill_matrix = true;
}
else
{
  $sparse /= 100;
}



if ($#ARGV < 1) {
  print "Usage: $0 n sparseness > file.matrix\n";
  print "       n is the integer size of the matrix\n";
  print "       sparseness is a percentage\n";
  exit -1; 
}


srand(time ^ $$ ^ unpack "%L*", `ps axww | gzip`);
for ($i = 1; $i <= $n; $i++)
{
  for ($j = 1; $j <= $n; $j++)
  {
    if ($i == $j) {
      print rand(100) . " ";
    } else {
      $nonzero = rand(1);
      if ($fill_matrix)
      {
        $num = rand(100);
        $sign = rand(1);
        if ($sign > 0.5){
          print "-";
        }
        print $num;
        print " ";
      }
      else
      {
        if ($nonzero > $sparse)
        {
          $num = rand(100);
          $sign = rand(1);
          if ($sign > 0.5){
            print "-";
          }
          print $num;
          print " ";
        }
        else
        {
          print "0 ";
        }
      }
    }
  }
  print "\n";
}




