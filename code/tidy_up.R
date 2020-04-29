
## tidy up the source code of image demo
x = file.path("modify_summary.R")

# to console
tidy_source(x)

# to a file
f = tempfile()
tidy_source(x, blank = TRUE, file = f, width.cutoff = 72, wrap = TRUE, indent = 2)

## check the original code here and see the difference
file.show(x)
file.show(f)
