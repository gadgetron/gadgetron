#!/usr/bin/perl

$gadgetron_home = "/home/ifedata/mrprogs/gadgetron";

my $executable = "$gadgetron_home/bin/gadgetron";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $timestring = sprintf "%4d%02d%02d_%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec;

print "Time string: $timestring\n";


$ENV{'GADGETRON_HOME'} = $gadgetron_home;
$ENV{'LD_LIBRARY_PATH'} = "/usr/local/lib:/usr/local/cuda/lib64:/usr/local/cula/lib64:" . $gadgetron_home . "/lib";

$exe_command = "killall gadgetron";
system($exe_command);

$exe_command = "mkdir -p log";
system($exe_command);

$exe_command = "nohup $executable > log/gadgetron_log_$timestring" . ".txt &" ;
system($exe_command);
