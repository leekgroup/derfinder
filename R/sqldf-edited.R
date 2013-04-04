# important function here is read.csv.sql
# arguments and return are exactly the same as read.csv.sql() in the sqldf package
# the only difference is that Tcl/Tk is NOT ALLOWED.

read.csv.sql <- function (file, sql = "select * from file", header = TRUE, sep = ",", 
    row.names, eol, skip, filter, nrows, field.types, comment.char, 
    dbname = tempfile(), drv = "SQLite", ...) 
{
	require(proto)
    file.format <- list(header = header, sep = sep)
    if (!missing(eol)) 
        file.format <- append(file.format, list(eol = eol))
    if (!missing(row.names)) 
        file.format <- append(file.format, list(row.names = row.names))
    if (!missing(skip)) 
        file.format <- append(file.format, list(skip = skip))
    if (!missing(filter)) 
        file.format <- append(file.format, list(filter = filter))
    if (!missing(nrows)) 
        file.format <- append(file.format, list(nrows = nrows))
    if (!missing(field.types)) 
        file.format <- append(file.format, list(field.types = field.types))
    if (!missing(comment.char)) 
        file.format <- append(file.format, list(comment.char = comment.char))
    pf <- parent.frame()
    if (missing(file) || is.null(file) || is.na(file)) 
        file <- ""
    tf <- NULL
    if (substring(file, 1, 7) == "http://" || substring(file, 
        1, 6) == "ftp://") {
        tf <- tempfile()
        on.exit(unlink(tf), add = TRUE)
        download.file(file, tf, mode = "wb")
        file <- tf
    }
    p <- proto(pf, file = file(file))
    p <- do.call(proto, list(pf, file = file(file)))
    sqldf(sql, envir = p, file.format = file.format, dbname = dbname, 
        drv = drv, ...)
}


####################################################################################
strapply <- function (X, pattern, FUN = function(x, ...) x, backref = NULL, 
    ..., empty = NULL, ignore.case = FALSE, perl = TRUE, engine = "R", 
    simplify = FALSE, USE.NAMES = FALSE, combine = c) 
{
    combine <- match.funfn(combine)
    stopifnot(!missing(pattern))
    pattern <- as.character(pattern)
    if (is.proto(FUN) || perl) 
        engine <- "R"
    if (identical(engine, "R")) 
        return(ostrapply(X = X, ignore.case = ignore.case, pattern = pattern, 
            FUN = FUN, backref = backref, ..., empty = empty, 
            perl = perl, simplify = simplify, USE.NAMES = USE.NAMES, 
            combine = combine))
    if (is.proto(FUN)) {
    }
    else if (is.character(FUN)) {
        FUN.orig <- FUN
        FUN <- function(...) FUN.orig
    }
    else if (is.list(FUN)) {
        values.replacement <- FUN
        names.replacement <- names(FUN)
        FUN <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    else {
        FUN <- match.funfn(FUN)
    }
    ff <- function(x) {
        s <- strapply1(x, pattern, backref, ignore.case)
        if (length(s) == 0 && !is.null(empty)) 
            s <- matrix(empty, 1)
        L <- lapply(seq_len(ncol(s)), function(j) {
            combine(do.call(FUN, as.list(s[, j])))
        })
        do.call("c", L)
    }
    result <- sapply(X, ff, simplify = is.logical(simplify) && 
        simplify, USE.NAMES = USE.NAMES)
    if (is.logical(simplify)) 
        result
    else do.call(match.funfn(simplify), result)
}

################################################################################################

match.funfn <- function (FUN, descend = TRUE) 
{
    if (is.function(FUN)) 
        return(FUN)
    if (inherits(FUN, "formula")) 
        return(as.function.formula(FUN))
    if (!(is.character(FUN) && length(FUN) == 1 || is.symbol(FUN))) {
        FUN <- eval.parent(substitute(substitute(FUN)))
        if (!is.symbol(FUN)) 
            stop(gettextf("'%s' is not a function, character or symbol", 
                deparse(FUN)), domain = NA)
    }
    envir <- parent.frame(2)
    if (descend) 
        FUN <- get(as.character(FUN), mode = "function", envir = envir)
    else {
        FUN <- get(as.character(FUN), mode = "any", envir = envir)
        if (!is.function(FUN)) 
            stop(gettextf("found non-function '%s'", FUN), domain = NA)
    }
    return(FUN)
}

################################################################################################

ostrapply <- function (X, pattern, FUN = function(x, ...) x, ignore.case = FALSE, 
    ..., empty = NULL, simplify = FALSE, USE.NAMES = FALSE, combine = c) 
{
    here <- environment()
    combine <- match.funfn(combine)
    if (is.character(FUN)) {
        FUN.orig <- FUN
        FUN <- function(...) FUN.orig
    }
    else if (is.list(FUN)) {
        values.replacement <- FUN
        names.replacement <- names(FUN)
        here$FUN <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    p <- if (is.proto(FUN)) {
        FUN$X <- X
        FUN$pattern <- pattern
        FUN$simplify <- simplify
        FUN$USE.NAMES <- USE.NAMES
        FUN$combine <- combine
        proto(pre = function(this) {
            this$first <- TRUE
            this$v <- NULL
            if (!is.null(FUN[["pre"]])) 
                FUN$pre()
        }, fun = function(this, ...) {
            FUN$count <- this$count
            this$v <- if (this$first) 
                combine(FUN$fun(...))
            else c(this$v, combine(FUN$fun(...)))
            this$first <- FALSE
        }, post = function(this) {
            if (this$first) 
                this$v <- NULL
            if (!is.null(FUN[["post"]])) 
                FUN$post()
        }, )
    }
    else {
        FUN <- match.funfn(FUN)
        proto(pre = function(this) {
            this$first <- TRUE
            this$v <- NULL
        }, fun = function(this, ...) {
            this$v <- if (this$first) 
                combine(FUN(...))
            else c(this$v, combine(FUN(...)))
            this$first <- FALSE
        }, post = function(this) {
            if (this$first) 
                this$v <- NULL
        })
    }
    ff <- function(x, ...) {
        gsubfn(pattern, p, x, ignore.case = ignore.case, 
            ...)
        p$v
    }
    result <- sapply(X, ff, ..., simplify = isTRUE(simplify), 
        USE.NAMES = USE.NAMES)
    if (is.logical(simplify)) 
        result
    else {
        do.call(match.funfn(simplify), result)
    }
}

################################################################################################

gsubfn <- function (pattern, replacement, x, backref, USE.NAMES = FALSE, 
    ignore.case = FALSE, env = parent.frame(), ...) 
{
	engine <- "R"
    R.engine <- identical(engine, "R")
    here <- environment()
    if (missing(replacement)) 
        here$replacement <- function(...) eval(parse(text = paste(..., 
            sep = "")), env)
    if (is.character(replacement)) {
        if (R.engine) 
            return(base::gsub(pattern, replacement, x, ...))
        else { stop("trying to use tcltk.  don't do that.")        }
    }
    if (is.list(replacement)) {
        values.replacement <- replacement
        names.replacement <- names(replacement)
        here$replacement <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    if (missing(pattern)) 
        pattern <- "[$]([[:alpha:]][[:alnum:].]*)|`([^`]+)`"
    pattern <- as.character(pattern)
    e <- NULL
    if (!inherits(replacement, "formula") && !is.function(replacement)) {
        e <- replacement
        e$pattern <- pattern
        e$x <- x
        e$backref <- if (!missing(backref)) 
            backref
        e$USE.NAMES <- USE.NAMES
        e$env <- env
        dots <- list(...)
        if (!is.null(names(dots))) {
            nam <- names(dots)
            for (n in nam[nam != ""]) assign(n, dots[[n]], e)
        }
        e$replacement <- function(this, ...) {
            this$count <- this$count + 1
            this$match <- c(...)
            this$fun(...)
        }
        here$replacement <- e$replacement
    }
    here$replacement <- match.funfn(replacement)
    if (missing(backref) || is.null(backref)) {
        noparen <- base::gsub("\\\\.", "", pattern)
        noparen <- base::gsub("\\[[^\\]]*\\]", "", noparen)
        backref <- -nchar(base::gsub("[^(]", "", noparen))
    }
    if (names(formals(here$replacement))[[1]] == "&") {
        backref <- abs(backref)
        if (!is.null(e)) 
            e$backref <- backref
    }
    j <- (identical(engine, "R") && !is.null(backref) && backref >= 
        0) + abs(backref)
    i <- if (!R.engine && backref >= 0) 
        0
    else 1
    j <- max(i, j)
    stopifnot(is.character(pattern), is.character(x), is.function(replacement))
    force(env)
    gsub.function <- function(x) {
        if (R.engine && !is.null(backref) && backref >= 0) {
            pattern <- paste("(", pattern, ")", sep = "")
        }
        if (!is.null(e)) {
            e$count <- 0
            if ("pre" %in% ls(e)) 
                e$pre()
        }
        repl <- function(i, j) {
            rs <- paste("\\", seq(i, j), collapse = "\002", sep = "")
            rs <- paste("\001\002", rs, "\001", sep = "")
            if (R.engine) 
                tryCatch(base::gsub(pattern, rs, x, ignore.case = ignore.case, 
                  ...), error = function(x) if (j > i) 
                  repl(i, j - 1)
                else stop(x))
            else { stop("trying to use tcltk. don't do that (2)")}
        }
        z <- repl(i, j)
        z <- strsplit(z, "\001")[[1]]
        f <- function(s) {
            if (nchar(s) > 0 && substring(s, 1, 1) == "\002") {
                s <- sub("\002$", "\002\002", s)
                L <- as.list(strsplit(s, "\002")[[1]][-1])
                do.call(replacement, L)
            }
            else s
        }
        z <- paste(sapply(z, f), collapse = "")
        if (!is.null(e) && "post" %in% ls(e)) 
            e$post()
        z
    }
    sapply(x, gsub.function, USE.NAMES = USE.NAMES)
}

################################################################################################

sqldf <- function (x, stringsAsFactors = FALSE, row.names = FALSE, envir = parent.frame(), 
    method = getOption("sqldf.method"), file.format = list(), 
    dbname, drv = getOption("sqldf.driver"), user, password = "", 
    host = "localhost", port, dll = getOption("sqldf.dll"), connection = getOption("sqldf.connection"), 
    verbose = isTRUE(getOption("sqldf.verbose"))) 
{
	require(DBI)
	#require(gsubfn)
	require(proto)
	require(chron)
	require(RSQLite)
    as.POSIXct.numeric <- function(x, ...) structure(x, class = c("POSIXct", 
        "POSIXt"))
    as.POSIXct.character <- function(x) structure(as.numeric(x), 
        class = c("POSIXct", "POSIXt"))
    as.Date.character <- function(x) structure(as.numeric(x), 
        class = "Date")
    as.Date2 <- function(x) UseMethod("as.Date2")
    as.Date2.character <- function(x) base:::as.Date.character(x)
    as.Date.numeric <- function(x, origin = "1970-01-01", ...) base::as.Date.numeric(x, 
        origin = origin, ...)
    as.dates.character <- function(x) structure(as.numeric(x), 
        class = c("dates", "times"))
    as.times.character <- function(x) structure(as.numeric(x), 
        class = "times")
    backquote.maybe <- function(nam) {
        if (drv == "h2") {
            nam
        }
        else if (drv == "mysql") {
            nam
        }
        else if (drv == "pgsql") {
            nam
        }
        else if (drv == "postgresql") {
            nam
        }
        else {
            if (regexpr(".", nam, fixed = TRUE)) {
                paste("`", nam, "`", sep = "")
            }
            else nam
        }
    }
    name__class <- function(data, ...) {
        if (is.null(data)) 
            return(data)
        cls <- sub(".*__([^_]+)|.*", "\\1", names(data))
        f <- function(i) {
            if (cls[i] == "") {
                data[[i]]
            }
            else {
                fun_name <- paste("as", cls[i], sep = ".")
                fun <- mget(fun_name, envir = environment(), 
                  mode = "function", ifnotfound = NA, inherits = TRUE)[[1]]
                if (identical(fun, NA)) 
                  data[[i]]
                else {
                  names(data)[i] <<- sub("__[^_]+$", "", names(data)[i])
                  fun(data[[i]])
                }
            }
        }
        data[] <- lapply(1:NCOL(data), f)
        data
    }
    colClass <- function(data, cls) {
        if (is.null(data)) 
            return(data)
        if (is.list(cls)) 
            cls <- unlist(cls)
        cls <- rep(cls, length = length(data))
        f <- function(i) {
            if (cls[i] == "") {
                data[[i]]
            }
            else {
                fun_name <- paste("as", cls[i], sep = ".")
                fun <- mget(fun_name, envir = environment(), 
                  mode = "function", ifnotfound = NA, inherits = TRUE)[[1]]
                if (identical(fun, NA)) 
                  data[[i]]
                else {
                  names(data)[i] <<- sub("__[^_]+$", "", names(data)[i])
                  fun(data[[i]])
                }
            }
        }
        data[] <- lapply(1:NCOL(data), f)
        data
    }
    overwrite <- FALSE
    request.open <- missing(x) && is.null(connection)
    request.close <- missing(x) && !is.null(connection)
    request.con <- !missing(x) && !is.null(connection)
    request.nocon <- !missing(x) && is.null(connection)
    dfnames <- fileobjs <- character(0)
    if (!is.list(method)) 
        method <- list(method, NULL)
    to.df <- method[[1]]
    to.db <- method[[2]]
    if (request.close || request.nocon) {
        on.exit({
            dbPreExists <- attr(connection, "dbPreExists")
            dbname <- attr(connection, "dbname")
            if (!missing(dbname) && !is.null(dbname) && dbname == 
                ":memory:") {
                if (verbose) {
                  cat("sqldf: dbDisconnect(connection)\n")
                }
                dbDisconnect(connection)
            } else if (!dbPreExists && drv == "sqlite") {
                if (verbose) {
                  cat("sqldf: dbDisconnect(connection)\n")
                  cat("sqldf: file.remove(dbname)\n")
                }
                dbDisconnect(connection)
                file.remove(dbname)
            } else {
                for (nam in dfnames) {
                  nam2 <- backquote.maybe(nam)
                  if (verbose) {
                    cat("sqldf: dbRemoveTable(connection, ", 
                      nam2, ")\n")
                  }
                  dbRemoveTable(connection, nam2)
                }
                for (fo in fileobjs) {
                  if (verbose) {
                    cat("sqldf: dbRemoveTable(connection, ", 
                      fo, ")\n")
                  }
                  dbRemoveTable(connection, fo)
                }
                if (verbose) {
                  cat("sqldf: dbDisconnect(connection)\n")
                }
                dbDisconnect(connection)
            }
        }, add = TRUE)
        if (request.close) {
            if (identical(connection, getOption("sqldf.connection"))) 
                options(sqldf.connection = NULL)
            return()
        }
    }
    if (request.open || request.nocon) {
        if (is.null(drv) || drv == "") {
            drv <- if ("package:RPostgreSQL" %in% search()) {
                "PostgreSQL"
            }
            else if ("package:RpgSQL" %in% search()) {
                "pgSQL"
            }
            else if ("package:RMySQL" %in% search()) {
                "MySQL"
            }
            else if ("package:RH2" %in% search()) {
                "H2"
            }
            else "SQLite"
        }
        drv <- sub("^[Rr]", "", drv)
        pkg <- paste("R", drv, sep = "")
        if (verbose) {
            if (!is.loaded(pkg)) 
                cat("sqldf: library(", pkg, ")\n", sep = "")
            library(pkg, character.only = TRUE)
        }
        else library(pkg, character.only = TRUE)
        drv <- tolower(drv)
        if (drv == "mysql") {
            if (verbose) 
                cat("sqldf: m <- dbDriver(\"MySQL\")\n")
            m <- dbDriver("MySQL")
            if (missing(dbname) || is.null(dbname)) {
                dbname <- getOption("RMySQL.dbname")
                if (is.null(dbname)) 
                  dbname <- "test"
            }
            connection <- if (missing(dbname) || dbname == ":memory:") {
                dbConnect(m)
            }
            else dbConnect(m, dbname = dbname)
            dbPreExists <- TRUE
        }
        else if (drv == "postgresql") {
            if (verbose) 
                cat("sqldf: m <- dbDriver(\"PostgreSQL\")\n")
            m <- dbDriver("PostgreSQL")
            if (missing(user) || is.null(user)) {
                user <- getOption("sqldf.RPostgreSQL.user")
                if (is.null(user)) 
                  user <- "postgres"
            }
            if (missing(password) || is.null(password)) {
                password <- getOption("sqldf.RPostgreSQL.password")
                if (is.null(password)) 
                  password <- "postgres"
            }
            if (missing(dbname) || is.null(dbname)) {
                dbname <- getOption("sqldf.RPostgreSQL.dbname")
                if (is.null(dbname)) 
                  dbname <- "test"
            }
            if (missing(host) || is.null(host)) {
                host <- getOption("sqldf.RPostgreSQL.host")
                if (is.null(host)) 
                  host <- "localhost"
            }
            if (missing(port) || is.null(port)) {
                port <- getOption("sqldf.RPostgreSQL.port")
                if (is.null(port)) 
                  port <- 5432
            }
            connection <- dbConnect(m, user = user, password, 
                dbname = dbname, host = host, port = port)
            if (verbose) 
                cat(sprintf("sqldf: connection <- dbConnect(m, user='%s', password=<...>, dbname = '%s', host = '%s', port = '%s')\n", 
                  user, dbname, host, port))
            dbPreExists <- TRUE
        }
        else if (drv == "pgsql") {
            if (verbose) 
                cat("sqldf: m <- dbDriver(\"pgSQL\")\n")
            m <- dbDriver("pgSQL")
            if (missing(dbname) || is.null(dbname)) {
                dbname <- getOption("RpgSQL.dbname")
                if (is.null(dbname)) 
                  dbname <- "test"
            }
            connection <- dbConnect(m, dbname = dbname)
            dbPreExists <- TRUE
        }
        else if (drv == "h2") {
            if (verbose) 
                cat("sqldf: m <- dbDriver(\"H2\")\n")
            m <- dbDriver("H2")
            if (missing(dbname) || is.null(dbname)) 
                dbname <- ":memory:"
            dbPreExists <- dbname != ":memory:" && file.exists(dbname)
            connection <- if (missing(dbname) || is.null(dbname) || 
                dbname == ":memory:") {
                dbConnect(m, "jdbc:h2:mem:", "sa", "")
            }
            else {
                jdbc.string <- paste("jdbc:h2", dbname, sep = ":")
                dbConnect(m, jdbc.string)
            }
        }
        else {
            if (verbose) 
                cat("sqldf: m <- dbDriver(\"SQLite\")\n")
            m <- dbDriver("SQLite")
            if (missing(dbname) || is.null(dbname)) 
                dbname <- ":memory:"
            dbPreExists <- dbname != ":memory:" && file.exists(dbname)
            dll <- getOption("sqldf.dll")
            if (length(dll) != 1 || identical(dll, FALSE) || 
                nchar(dll) == 0) {
                dll <- FALSE
            }
            else {
                if (dll == basename(dll)) 
                  dll <- Sys.which(dll)
            }
            options(sqldf.dll = dll)
            if (!identical(dll, FALSE)) {
                if (verbose) {
                  cat("sqldf: connection <- dbConnect(m, dbname = \"", 
                    dbname, "\", loadable.extensions = TRUE\n", 
                    sep = "")
                  cat("sqldf: library(RSQLite.extfuns)\n")
                  cat("sqldf: select load_extension('", dll, 
                    "')\n", sep = "")
                }
                connection <- dbConnect(m, dbname = dbname, loadable.extensions = TRUE)
                library("RSQLite.extfuns", character.only = TRUE)
                s <- sprintf("select load_extension('%s')", dll)
                dbGetQuery(connection, s)
            }
            else {
                if (verbose) {
                  cat("sqldf: connection <- dbConnect(m, dbname = \"", 
                    dbname, "\")\n", sep = "")
                }
                connection <- dbConnect(m, dbname = dbname)
            }
            if (verbose) 
                cat("sqldf: init_extensions(connection)\n")
            library(RSQLite.extfuns)
            init_extensions(connection)
        }
        attr(connection, "dbPreExists") <- dbPreExists
        if (missing(dbname) && drv == "sqlite") 
            dbname <- ":memory:"
        attr(connection, "dbname") <- dbname
        if (request.open) {
            options(sqldf.connection = connection)
            return(connection)
        }
    }
    if (request.con) {
        drv <- if (inherits(connection, "PostgreSQLConnection")) 
            "PostgreSQL"
        else if (inherits(connection, "pgSQLConnection")) 
            "pgSQL"
        else if (inherits(connection, "MySQLConnection")) 
            "MySQL"
        else if (inherits(connection, "H2Connection")) 
            "H2"
        else "SQLite"
        drv <- tolower(drv)
        dbPreExists <- attr(connection, "dbPreExists")
    }
    engine <- "R" ## removing option to require tcltk
    words. <- words <- if (engine == "tcl") {
        strapplyc(x, "[[:alnum:]._]+")
    }
    else strapply(x, "[[:alnum:]._]+", engine = "R")
    if (length(words) > 0) 
        words <- unique(unlist(words))
    is.special <- sapply(mget(words, envir, "any", NA, inherits = TRUE), 
        function(x) is.data.frame(x) + 2 * inherits(x, "file"))
    dfnames <- words[is.special == 1]
    for (i in seq_along(dfnames)) {
        nam <- dfnames[i]
        if (dbPreExists && !overwrite && dbExistsTable(connection, 
            nam)) {
            dfnames <- head(dfnames, i - 1)
            stop(paste("sqldf:", "table", nam, "already in", 
                dbname, "\n"))
        }
        DF <- as.data.frame(get(nam, envir))
        if (!is.null(to.db) && is.function(to.db)) 
            DF <- to.db(DF)
        nam2 <- backquote.maybe(nam)
        if (verbose) 
            cat("sqldf: dbWriteTable(connection, '", nam2, "', ", 
                nam, ", row.names = ", row.names, ")\n", sep = "")
        dbWriteTable(connection, nam2, DF, row.names = row.names)
    }
    fileobjs <- if (is.null(file.format)) {
        character(0)
    }
    else {
        eol <- if (.Platform$OS == "windows") 
            "\r\n"
        else "\n"
        words[is.special == 2]
    }
    for (i in seq_along(fileobjs)) {
        fo <- fileobjs[i]
        Filename <- summary(get(fo, envir))$description
        if (dbPreExists && !overwrite && dbExistsTable(connection, 
            Filename)) {
            fileobjs <- head(fileobjs, i - 1)
            stop(paste("sqldf:", "table", fo, "from file", Filename, 
                "already in", dbname, "\n"))
        }
        args <- c(list(conn = connection, name = fo, value = Filename), 
            modifyList(list(eol = eol, comment.char = "", quote = "\""), 
                file.format))
        args <- modifyList(args, as.list(attr(get(fo, envir), 
            "file.format")))
        filter <- args$filter
        if (!is.null(filter)) {
            args$filter <- NULL
            Filename.tmp <- tempfile()
            args$value <- Filename.tmp
            filter.subs <- filter[-1]
            if (length(filter.subs) > 0) {
                filter.subs <- filter.subs[sapply(names(filter.subs), 
                  nzchar)]
            }
            filter.nms <- names(filter.subs)
            filter.tempfiles <- sapply(filter.nms, tempfile)
            cmd <- filter[[1]]
            for (nm in filter.nms) {
                cat(filter.subs[[nm]], file = filter.tempfiles[[nm]])
                cmd <- gsub(nm, filter.tempfiles[[nm]], cmd, 
                  fixed = TRUE)
            }
            cmd <- if (nchar(Filename) > 0) 
                sprintf("%s < \"%s\" > \"%s\"", cmd, Filename, 
                  Filename.tmp)
            else sprintf("%s > \"%s\"", cmd, Filename.tmp)
            # if (.Platform$OS == "windows") {
                # cmd <- paste("cmd /c", cmd)
                # if (FALSE) {
                  # key <- "SOFTWARE\\R-core"
                  # show.error.messages <- getOption("show.error.message")
                  # options(show.error.messages = FALSE)
                  # reg <- try(readRegistry(key, maxdepth = 3)$Rtools$InstallPath)
                  # reg <- NULL
                  # options(show.error.messages = show.error.messages)
                  # if (!is.null(reg) && !inherits(reg, "try-error")) {
                    # Rtools.path <- file.path(reg, "bin", fsep = "\\")
                    # path <- Sys.getenv("PATH")
                    # on.exit(Sys.setenv(PATH = path), add = TRUE)
                    # path.new <- paste(path, Rtools.path, sep = ";")
                    # Sys.setenv(PATH = path.new)
                  # }
                # }
            # }
            ## causes problems in compiling the R package - this isn't meant for windows anyway.
            if (verbose) 
                cat("sqldf: system(\"", cmd, "\")\n", sep = "")
            system(cmd)
            for (fn in filter.tempfiles) file.remove(fn)
        }
        if (verbose) 
            cat("sqldf: dbWriteTable(", toString(args), ")\n")
        do.call("dbWriteTable", args)
    }
    if (drv == "sqlite" || drv == "mysql" || drv == "postgresql") {
        for (xi in x) {
            if (verbose) {
                cat("sqldf: dbGetQuery(connection, '", xi, "')\n", 
                  sep = "")
            }
            rs <- dbGetQuery(connection, xi)
        }
    }
    else {
        for (i in seq_along(x)) {
            if (length(words.[[i]]) > 0) {
                dbGetQueryWords <- c("select", "show", "call", 
                  "explain", "with")
                if (tolower(words.[[i]][1]) %in% dbGetQueryWords) {
                  if (verbose) {
                    cat("sqldf: dbGetQuery(connection, '", x[i], 
                      "')\n", sep = "")
                  }
                  rs <- dbGetQuery(connection, x[i])
                }
                else {
                  if (verbose) {
                    cat("sqldf: dbSendUpdate:", x[i], "\n")
                  }
                  rs <- dbSendUpdate(connection, x[i])
                }
            }
        }
    }
    if (is.null(to.df)) 
        to.df <- "auto"
    if (is.function(to.df)) 
        return(to.df(rs))
    if (identical(to.df, "raw")) 
        return(rs)
    if (identical(to.df, "name__class")) 
        return(do.call("name__class", list(rs)))
    if (!identical(to.df, "nofactor") && !identical(to.df, "auto")) {
        return(do.call("colClass", list(rs, to.df)))
    }
    row_names_name <- grep("row[_.]names", names(rs), value = TRUE)
    if (length(row_names_name) > 1) 
        warning(paste("ambiguity regarding row names:", row_names_name))
    row_names_name <- row_names_name[1]
    rs <- if (!is.na(row_names_name)) {
        if (identical(row.names, FALSE)) {
            rs[names(rs) != row_names_name]
        }
        else {
            rn <- rs[[row_names_name]]
            rs <- rs[names(rs) != row_names_name]
            if (all(regexpr("^[[:digit:]]*$", rn) > 0)) 
                rn <- as.integer(rn)
            rownames(rs) <- rn
            rs
        }
    }
    else rs
    tab <- do.call("rbind", lapply(dfnames, function(dfname) {
        df <- get(dfname, envir)
        nms <- names(df)
        do.call("rbind", lapply(seq_along(df), function(j) {
            column <- df[[j]]
            cbind(dfname, nms[j], toString(class(column)), toString(levels(column)))
        }))
    }))
    tabu <- unique(tab[, -1, drop = FALSE])
    dup <- unname(tabu[duplicated(tabu[, 1]), 1])
    auto <- function(i) {
        cn <- colnames(rs)[[i]]
        if (!cn %in% dup && (ix <- match(cn, tab[, 2], nomatch = 0)) > 
            0) {
            df <- get(tab[ix, 1], envir)
            if (inherits(df[[cn]], "ordered")) {
                if (identical(to.df, "auto")) {
                  u <- unique(rs[[i]])
                  levs <- levels(df[[cn]])
                  if (all(u %in% levs)) 
                    return(factor(rs[[i]], levels = levels(df[[cn]]), 
                      ordered = TRUE))
                  else return(rs[[i]])
                }
                else return(rs[[i]])
            }
            else if (inherits(df[[cn]], "factor")) {
                if (identical(to.df, "auto")) {
                  u <- unique(rs[[i]])
                  levs <- levels(df[[cn]])
                  if (all(u %in% levs)) 
                    return(factor(rs[[i]], levels = levels(df[[cn]])))
                  else return(rs[[i]])
                }
                else return(rs[[i]])
            }
            else if (inherits(df[[cn]], "POSIXct")) 
                return(as.POSIXct(rs[[i]]))
            else if (inherits(df[[cn]], "times")) 
                return(as.times.character(rs[[i]]))
            else {
                asfn <- paste("as", class(df[[cn]]), sep = ".")
                asfn <- match.fun(asfn)
                return(asfn(rs[[i]]))
            }
        }
        if (stringsAsFactors && is.character(rs[[i]])) 
            factor(rs[[i]])
        else rs[[i]]
    }
    rs2 <- lapply(seq_along(rs), auto)
    rs[] <- rs2
    rs
}

