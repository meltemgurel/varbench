library(data.table)
poss <- function(v){
      m <- c()
      p <- v[1]
      while(is.finite(p)){
        m <- append(m, p)
        p <- min(v[v >= p + 101])
      }
      return(m)
    }
