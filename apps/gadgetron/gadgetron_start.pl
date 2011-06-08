#!/usr/bin/perl

$home_dir = $ENV{"HOME"};
$gadgetron_home = $home_dir . "/mrprogs/gadgetron";

my $executable = "$gadgetron_home/bin/gadgetron";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $timestring = sprintf "%4d%02d%02d_%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec;

print "Time string: $timestring\n";


$ENV{'GADGETRON_HOME'} = $gadgetron_home;
$ENV{'LD_LIBRARY_PATH'} = "/usr/local/lib:/usr/local/cuda/lib64:/usr/local/cula/lib64:" . $gadgetron_home . "/lib";

$exe_command = "killall -9 gadgetron";
system($exe_command);
sleep(1);

$exe_command = "mkdir -p log";
system($exe_command);

$logfilename = "log/gadgetron_log_$timestring" . ".txt";
$exe_command = "nohup $executable > $logfilename 2> $logfilename < /dev/null &" ;
system($exe_command);

sleep(1);

