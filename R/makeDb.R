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





#'Create SQLite database from text file
#'
#'Dumps the contents of a table (saved as a text file) into a SQLite database,
#'performing some filtering along the way.
#'
#'
#'@param dbfile Character string giving the file name/location of the database
#'to be created.  Generally ends in \code{.db}.
#'@param textfile The text file containing the table to be dumped into
#'\code{dbfile}.
#'@param tablename Character string containing name to give the table inside
#'\code{dbfile}.
#'@param sep The separator used in \code{textfile}.  The tornado pipeline
#'creates tab-separated text files, so \code{"\t"} is the default.
#'@param cutoff Rows in \code{textfile} must have at least one entry (not
#'counting the first column, which is assumed to hold genomic position) greater
#'than \code{cutoff} to be included in \code{dbfile}.
#'@return No return, but writes the file \code{dbfile} containing table
#'\code{tablename} by filtering \code{textfile} according to \code{cutoff}.
#'@note The workhorse of this function is a modified version of
#'\code{\link{read.csv.sql}}, found in the \code{sqldf} package.
#'@author Alyssa Frazee
#'@export

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
