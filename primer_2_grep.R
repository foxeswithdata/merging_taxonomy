primer_to_grep <- function(primer, forward =TRUE){
  
  primer <- stringr::str_replace(primer, "N", "[ACGT]")
  primer <- stringr::str_replace(primer, "R", "[AG]")
  primer <- stringr::str_replace(primer, "Y", "[CT]")
  primer <- stringr::str_replace(primer, "M", "[AC]")
  primer <- stringr::str_replace(primer, "K", "[GT]")
  primer <- stringr::str_replace(primer, "S", "[CG]")
  primer <- stringr::str_replace(primer, "W", "[AT]")
  primer <- stringr::str_replace(primer, "B", "[CGT]")
  primer <- stringr::str_replace(primer, "D", "[AGT]")
  primer <- stringr::str_replace(primer, "H", "[ACT]")
  primer <- stringr::str_replace(primer, "V", "[ACG]")
  
  return(primer)
}
