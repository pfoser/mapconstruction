function ts = timestamp()
% 
% TIMESTAMP returns the current system time.
% The output is type string and without EOL.
% (LB 9'02)

TempTime=clock;
ts = ['(' num2str(TempTime(4),'%02.0f') ':' num2str(TempTime(5),'%02.0f') ':' num2str(TempTime(6),'%02.0f') ')'];
