primer_to_grep <- function(primer, forward =TRUE){
  
  primer <- stringr::str_replace_all(primer, "N", "[ACGT]")
  primer <- stringr::str_replace_all(primer, "R", "[AG]")
  primer <- stringr::str_replace_all(primer, "Y", "[CT]")
  primer <- stringr::str_replace_all(primer, "M", "[AC]")
  primer <- stringr::str_replace_all(primer, "K", "[GT]")
  primer <- stringr::str_replace_all(primer, "S", "[CG]")
  primer <- stringr::str_replace_all(primer, "W", "[AT]")
  primer <- stringr::str_replace_all(primer, "B", "[CGT]")
  primer <- stringr::str_replace_all(primer, "D", "[AGT]")
  primer <- stringr::str_replace_all(primer, "H", "[ACT]")
  primer <- stringr::str_replace_all(primer, "V", "[ACG]")
  
  return(primer)
}
