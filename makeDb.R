## makeDb():
## arguments:
## --dbfile: character string giving the name/location of the database file you want to create
## --textfile: character string giving the name/location of the tab-separated coverage table (output of preprocessing pipeline)
## --tablename: character string giving the name of the table that will be created within dbfile
## --sep: the separator used in textfile.  (our software gives tab-separated coverage files).
## --cutoff: only rows of textfile with at least one sample having coverage bigger than cutoff will be dumped into dbfile.
## return:
## --no return, but creates an sql database dbfile containing table tablename from textfile, including only rows meeting the condition specifed by cutoff
## --the database created can be used by genominator so you don't have to load a giant table into memory.



makeDb <- function(dbfile, textfile, tablename, sep = "\t", cutoff = 5){

	cat(file=dbfile) #create empty file
	column.names = as.character(as.matrix(read.table(textfile,sep=sep,nrows=1,header=F))) #get column names
	#print(column.names)
	# create the sql statement to be used in creating database:
	for(i in 2:length(column.names)){ 
	  if(i==2) where.statement = paste(gsub("\\.","_",column.names[i]),">",cutoff,"OR")
	  else if(i==length(column.names)) where.statement = paste(where.statement,
	  	gsub("\\.","_",column.names[i]),">",cutoff)
  	  else where.statement = paste(where.statement,gsub("\\.","_",column.names[i]),">",cutoff,"OR")
  	  } #(note that "." is not an acceptable character in sql column names - is automatically replaced with "_" - so we do the same in our sql statement)
	tablename.statement = paste("main",tablename,sep=".")
	sql.statement = paste("create table",tablename.statement,"as select * from file where",where.statement)
	#print(where.statement)
	#print(sql.statement)
	read.csv.sql(textfile, sql=sql.statement, dbname=dbfile,sep=sep) #this is where all the action is - this creates the database, this is the function in our modified sqldf script.
	print(paste("Wrote database file",dbfile,"containing table",tablename)) # so that users know this function did something and created a file in their system.
}
