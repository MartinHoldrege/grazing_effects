# paths (put here, and sourced elsewhere so paths in each script don't
# need to be updated)

# path to where large files are stored, that aren't in working directory

if(dir.exists("E:/USGS")) {
  path_sw <- "E:/USGS/large_files/stepwat"
  path_large <- "E:/USGS/large_files"
} else  if(dir.exists("D:/USGS")) {
    path_sw <- "D:/USGS/large_files/stepwat"
    path_large <- "D:/USGS/large_files"
} else {
  warning('directory not found')
}

