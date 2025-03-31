packages <- c('tidyverse', 'dunn.test', 'rcompanion')

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

invisible(lapply(packages, install_if_missing))

lapply(packages, library, character.only = TRUE)

options(timeout = max(3600, getOption('timeout')))
names <- c(la_esperanza = 'La Esperanza',rio_claro = 'Rio Claro',
           T1 = 'T1',T2 = 'T2',T3 = 'T3',T4 = 'T4',T0 = 'T0',
           '1' = '1','2' = '2','3' = '3', 
           '2022-2023' = '2022-2023', '2023-2024' = '2023-2024')

cld <- function(variable, grupo) {
  
  dunn_result <- dunn.test(variable, grupo, method = 'bonferroni')
  
  if (length(unique(grupo)) > 2) {
    cld_c <- cldList(comparison = dunn_result$comparisons, 
                     p.value = dunn_result$P.adjusted, 
                     threshold = .05/2) |>
      rename(grupo = Group,
             cld = Letter) |>
      select(-MonoLetter) |>
      mutate(grupo = ifelse(grupo == 'T','T0',grupo))
    
  } else {
    
    dunn_result <- dunn.test(variable, grupo, method = 'bonferroni')
    
    if (dunn_result$P.adjusted <= 0.05/2) {cld_l = c('a','b')} else {cld_l = c('a','a')}
    
    cld_c <- tibble(grupo = unique(grupo), cld = cld_l)
    
  }
  
  return(cld_c)
}
