%Author: KB
%Purpose: For giving the date and time of execution of a script. Insert to 
%keep track of multiple executions of a script on the same machine or
%across several machines.

function runtimeTimeStamp()
dateAndTime=datetime('now','TimeZone','local','Format','dd-MMM-y HH:mm:ss Z');
disp(dateAndTime);
end