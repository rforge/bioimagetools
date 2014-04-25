.status<-function(status)
{
if (is.null(status))
{
status=0
cat("-")
return(0)
}
switch(status,
"0" = cat("\b-"),
"1" = cat("\b\\"),
"2" = cat("\b|"),
"3" = cat("\b/"),
)
status=status+1
if (status==4)status=0
return(status)
}
