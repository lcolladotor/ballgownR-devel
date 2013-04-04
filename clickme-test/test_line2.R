library(clickme)
library(colorspace)

### Custom line

## set path
set_root_path("/Users/lcollado/enigma2/ballgownR-devel/clickme-test/")



mammals_path <- system.file(file.path("examples", "one_zoom", "data", "mammals.tree"), package="clickme")
cities <- c("Boston", "NYC", "Philadelphia")
n <- 30
df2 <- data.frame(line = rep(cities, each = n), x = rep(1:n, length(cities)), y = c(sort(rnorm(n)), -sort(rnorm(n)),sort(rnorm(n))))

## Initial manual color test
col <- '["#393b79","#5254a3","#6b6ecf"]'
clickme(df2, "line_with_focus2", params=list(color=col), html="test.html")

## Test with groups
test <- data.frame(line=rep(letters[1:9], each=30), x=rep(1:30, 9), y=rnorm(30*9))
col2 <- paste0('["', paste0(rep(rainbow_hcl(3, start = 30, end = 300), each=3), collapse='","'), '"]')
clickme(test, "line_with_focus2", params=list(color=col2), html="test2.html")

## Testing whether you can show an "exon"
col <- '["#393b79","#5254a3","#6b6ecf","#6b6ecf"]'
df3 <- rbind(df2, data.frame(line=rep("exon", 5), x=11:15, y=-5))
clickme(df3, "line_with_focus2", params=list(color=col), html="test3.html")


## without specifying color param
clickme(df3, "line_with_focus2", html="test4.html")