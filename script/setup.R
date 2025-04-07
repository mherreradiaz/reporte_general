packages <- c('tidyverse', 'dunn.test', 'rcompanion','car')

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

cld_dunn <- function(variable, grupo) {
  
  dunn_result <- dunn.test(variable, grupo, method = 'bonferroni')
  
  if (length(unique(grupo)) > 2) {
    cld_c <- cldList(
      comparison = dunn_result$comparisons, 
      p.value = dunn_result$P.adjusted, 
      threshold = 0.05
    ) |>
      rename(grupo = Group,
                    cld = Letter) |>
      select(-MonoLetter)
    
  } else {
    cld_l <- ifelse(dunn_result$P.adjusted <= 0.05, c('a','b'), c('a','a'))
    cld_c <- tibble::tibble(grupo = unique(grupo), cld = cld_l)
  }
  
  cld_c <- cld_c |>
    mutate(grupo = ifelse(grupo == 'T', 'T0', as.character(grupo)))
  
  return(cld_c)
}

is.normal <- function(x, alpha = 0.05, detect_outliers = TRUE) {
  x <- na.omit(x)
  n <- length(x)

  if (n < 5) {
    warning(
      "El tamaño muestral es demasiado pequeño (n < 5). ",
      "Los tests de normalidad no son confiables en este rango."
    )
    return(NA)
  }

  has_outliers <- FALSE
  if (detect_outliers && n >= 10) {
    Q <- quantile(x, probs = c(0.25, 0.75))
    iqr <- IQR(x)
    lower <- Q[1] - 1.5 * iqr
    upper <- Q[2] + 1.5 * iqr
    has_outliers <- any(x < lower | x > upper)
  }

  test_name <- if (n <= 50) "Shapiro-Wilk" else 
    if (has_outliers) "Anderson-Darling" else "Kolmogorov-Smirnov"

  p_value <- switch(
    test_name,
    "Shapiro-Wilk" = shapiro.test(x)$p.value,
    "Anderson-Darling" = {
      if (!require(nortest, quietly = TRUE)) install.packages("nortest")
      nortest::ad.test(x)$p.value
    },
    "Kolmogorov-Smirnov" = ks.test(x, "pnorm", mean(x), sd(x))$p.value
  )

  is_normal <- if (n < 5) NA else p_value > alpha

  if (n >= 5) {
    message(
      sprintf(
        "Test: %s | p = %.4f | Outliers: %s | Conclusión: %snormal (n = %d)",
        test_name, 
        p_value,
        ifelse(has_outliers, "Sí", "No"),
        ifelse(is_normal, "", "NO "),
        n
      )
    )
  }
  
  return(is_normal)
}
is_all_true <- function(x) {
  !anyNA(x) && all(x)
}

dif_summary <- \(data,var,group,sup_group) {
  
  homocedasticidad <- data |> 
    group_by(!!!syms(sup_group)) |> 
    reframe(homocedasticidad = {
      x <- !!!syms(var)
      g <- as.factor(!!!syms(group))
      f_value <- leveneTest(x ~ g, data = pick(everything()))$"Pr(>F)"[1]
      f_value > 0.05
    })
  
  normalidad <- data |> 
    group_by(!!!syms(c(sup_group,group))) |> 
    reframe(normalidad = is.normal(!!!syms(var))) |> 
    suppressWarnings() |> 
    group_by(!!!syms(sup_group)) |> 
    reframe(normalidad = is_all_true(normalidad))
  
  test <- left_join(normalidad,homocedasticidad) |>
    mutate(test = case_when(normalidad == T & homocedasticidad == T ~ 'anova + tukey'
                            ,.default = 'kruskall-wallis + dunn'))
  
  data_test <- data |> 
    left_join(test)
  
  available_test <- data_test |> 
    pull(test) |> 
    unique()
  
  if ('kruskall-wallis + dunn' %in% available_test) {
    kw <- data_test |> 
      filter(test == 'kruskall-wallis + dunn') |> 
      group_by(!!!syms(sup_group)) |> 
      reframe({
        x = !!!syms(var)
        g = !!!syms(group)
        p_value <- kruskal.test(x ~ g, data = pick(everything()))$p.value
        tibble(p_value = p_value, significativo = p_value < 0.05)
      })
  } else kw <- data_test |> distinct(!!!syms(sup_group))
  
  if ('anova + tukey'  %in% available_test) {
    anova <- data_test |> 
      filter(test == 'anova + tukey') |> 
      group_by(!!!syms(sup_group)) |> 
      reframe(diferencias = {
        x = !!!syms(var)
        g = !!!syms(group)
        p_value <- summary(aov(x ~ g, data = pick(everything())))[[1]]$'Pr(>F)'[1]
        p_value < 0.05
      })
  } else anova <- data_test |> distinct(!!!syms(sup_group))
  
  diferencias <- left_join(kw,anova) |> 
    suppressMessages()
  
  data_diferencias <- data_test |>
    left_join(diferencias)
  
  available_test <- data_diferencias |> 
    filter(significativo == T) |> 
    pull(test) |> 
    unique()
  
  if (any(pull(diferencias,significativo))) {
    
    if ('kruskall-wallis + dunn' %in% available_test) {
      dunn <- data_diferencias |> 
        filter(significativo == T,
               test == 'kruskall-wallis + dunn') |> 
        group_by(!!!syms(sup_group)) |> 
        reframe(cld = cld_dunn(!!!syms(var),!!!syms(group))$cld,
                !!sym(group) := cld_dunn(!!!syms(var),!!!syms(group))$grupo) |> 
        select(all_of(c(sup_group,group)),cld)
    } else dunn <- data_diferencias |> distinct(!!!syms(c(sup_group,group)))
    
    if ('anova + tukey'  %in% available_test) {
      tukey <- data_diferencias |> 
        filter(significativo == T,
               test == 'anova + tukey') |> 
        group_by(!!!syms(sup_group)) |> 
        reframe({
          modelo <- aov(!!!syms(var) ~ !!!syms(group), data = pick(everything()))
          tukey <- TukeyHSD(modelo)
          cld_tukey <- (multcompLetters4(modelo, tukey)$tratamiento)$monospacedLetters
          tibble(!!!syms(group) := names(cld_tukey),cld = as.character(cld_tukey))
        })
    } else tukey <- data_diferencias |> distinct(!!!syms(c(sup_group,group)))
    
    grupos <- left_join(tukey,dunn) |> 
      suppressMessages() |> 
      mutate(cld = ifelse(is.na(cld),'a',cld))
    
    data_diferencias |> 
      left_join(grupos) |> 
      select(all_of(c(sup_group,group)),,normalidad,homocedasticidad,test,p_value,significativo,cld) |> 
      distinct()
    
  } else {
    
    data |> 
      select(all_of(c(sup_group,group))) |> 
      distinct() |> 
      mutate(cld = 'a')
  }
}

dif_summary_two <- \(data,var,group,sup_group) {
  
  normalidad <- data |> 
    group_by(!!!syms(c(sup_group,group))) |> 
    reframe(normalidad = is.normal(!!!syms(var))) |> 
    suppressWarnings() |> 
    group_by(!!!syms(sup_group)) |> 
    reframe(normalidad = is_all_true(normalidad))
  
  test <- normalidad |>
    mutate(test = case_when(normalidad == T ~ 'anova + tukey'
                            ,.default = 'kruskall-wallis + dunn'))
  
  data_test <- data |> 
    left_join(test)
  
  available_test <- data_test |> 
    pull(test) |> 
    unique()
  
  if ('kruskall-wallis + dunn' %in% available_test) {
    kw <- data_test |> 
      filter(test == 'kruskall-wallis + dunn') |> 
      group_by(!!!syms(sup_group),test) |> 
      reframe({
        x = !!!syms(var)
        g = !!!syms(group)
        p_value <- kruskal.test(x ~ g, data = pick(everything()))$p.value
        tibble(p_value = p_value, significativo = p_value < 0.05)
      })
  } else kw <- data_test |> distinct(!!!syms(sup_group))
  
  if ('anova + tukey'  %in% available_test) {
    anova <- data_test |> 
      filter(test == 'anova + tukey') |> 
      group_by(!!!syms(sup_group),test) |> 
      reframe(diferencias = {
        x = !!!syms(var)
        g = !!!syms(group)
        p_value <- summary(aov(x ~ g, data = pick(everything())))[[1]]$'Pr(>F)'[1]
        p_value < 0.05
      })
  } else anova <- data_test |> distinct(!!!syms(sup_group))
  
  diferencias <- left_join(kw,anova) |> 
    suppressMessages()
  
  data_test |> 
    distinct(!!!syms(c(sup_group,group))) |> 
    left_join(diferencias) |> 
    mutate(cld = case_when(significativo & temporada == '2022-2023'~'a',
                           significativo & temporada == '2023-2024'~'b',
                           .default = 'a'))
}

max_label <- \(var) {
  max_out <- max(boxplot.stats(var)$out)
  max_out <- ifelse(is.infinite(max_out),NA,max_out)
  
  max_stats <- max(boxplot.stats(var)$stats)
  
  return(max(c(max_out,max_stats),na.rm=T))
}

graficar_parametros <- \(dir,var,huerto,label) {
  
  cosecha_la_esperanza <- tibble(sitio = 'la_esperanza', fecha = as.Date('2022-12-12'), temporada = c('2022-2023','2023-2024'))
  cosecha_rio_claro <- tibble(sitio = 'rio_claro', fecha = c(as.Date('2022-12-21'),as.Date('2023-01-12')), temporada = c('2022-2023','2023-2024'))
  
  cosecha <- bind_rows(cosecha_la_esperanza,cosecha_rio_claro) |> 
    crossing(comparacion = factor(c('T1', 'T2', 'T3', 'T4'), levels = c('T1', 'T2', 'T3', 'T4')))
  
  fechas_mensuales <- seq.Date(from = as.Date("2022-09-01"), 
                               to   = as.Date("2023-05-01"), 
                               by   = "month")
  
  data <- read_rds(dir) |> 
    group_by(sitio,temporada,tratamiento,fecha) |> 
    reframe(mean_value = mean(!!sym(var),na.rm=T)) |> 
    mutate(tratamiento = factor(tratamiento,levels=paste0('T',0:4)))
  
  data_plot <- data |> 
    filter(tratamiento != 'T0') |> 
    mutate(comparacion = tratamiento) |> 
    bind_rows(
      data |> 
        filter(tratamiento == 'T0') |> 
        crossing(comparacion = factor(c('T1', 'T2', 'T3', 'T4'), levels = c('T1','T2','T3','T4')))
    ) |> 
    mutate(fecha = case_when(temporada == '2023-2024' ~ fecha - years(1),
                             .default = fecha))
  
  data_plot |> 
    filter(sitio == huerto) |> 
    ggplot(aes(fecha,mean_value, linetype = tratamiento, color = tratamiento)) +
    geom_line(linewidth = .8, alpha= .7) +
    geom_vline(xintercept = fechas_mensuales, linetype = "dotted", color = "grey30") +
    geom_vline(data = cosecha |> filter(sitio == huerto), aes(xintercept = fecha), 
               linewidth = .8, linetype = 'dashed', color = 'green4', alpha = .7) +
    scale_linetype_manual(values = c('T0' = 'dotted', 'T1' = 'solid', 'T2' = 'solid', 'T3' = 'solid', 'T4' = 'solid')) +
    scale_color_manual(values = c('T0' ="grey10",'T1' = "royalblue4",'T2' = "lightskyblue", 'T3' = "rosybrown2", 'T4' = "red3")) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(.1,0)) +
    labs(x = NULL, y = label, color = 'Tratamiento') +
    facet_grid(comparacion~temporada, scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_rect(fill = 'white'),
          legend.position = 'none',
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(linetype = "dashed"))
}

