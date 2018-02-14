BEGIN{
  found = 0;
}

/^subroutine/{
  found = 1;
  split($2, string, "(");
  name = string[1]".f90"
  system("rm -f "name);
  print "!######################################################################" >> name
  print $0 >> name
}

/integer function/ || /real\(8\) function/ || /complex\(8\) function/{
  found = 1;
  split($3, string, "(");
  name = string[1]".f90"
  system("rm -f "name);
  print "!######################################################################" >> name
  print $0 >> name
}

! /subroutine/ && ! / function / && ! /function$/{
  if (found == 1) {
    print $0 >> name
  }
}

/^end/{
  found = 0;
  print $0 >> name
  print "!######################################################################" >> name
}
