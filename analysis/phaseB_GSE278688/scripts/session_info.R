# session_info.R -- exact package versions for the frozen run.
cat("=== R ===\n")
cat(R.version.string, "\n")
cat("platform :", R.version$platform, "\n")
cat("nickname :", R.version$nickname, "\n\n")
cat("=== sessionInfo() ===\n")
print(sessionInfo())
cat("\n=== explicitly used packages ===\n")
for (p in c("data.table","Matrix","boot","statmod","jsonlite")) {
  v <- tryCatch(as.character(packageVersion(p)), error=function(e) "NOT INSTALLED")
  cat(sprintf("  %-12s %s\n", p, v))
}
