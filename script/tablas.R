source('script/setup.R')

data_apariencia <- read_rds('data/apariencia.rds')

data_produccion <- read_rds('data/produccion.rds')

data_brix <- read_rds('data/brix.rds')

data_daño <- read_rds('data/daño.rds') |> 
  mutate(daño = 100-(d_nulo/n)*100)

peso_total <- data_produccion |> 
  group_by(sitio,temporada,tratamiento) |> 
  reframe(mediana = median(peso_total,na.rm=T)) |> 
  group_by(sitio,temporada) |> 
  mutate(dif_T0 = round(mediana-first(mediana),2))

densidad <- data_produccion |> 
  group_by(sitio,temporada,tratamiento) |> 
  reframe(mediana = median(densidad,na.rm=T)) |> 
  group_by(sitio,temporada) |> 
  mutate(dif_T0 = round(mediana-first(mediana),2))

peso <- data_apariencia |> 
  group_by(sitio,temporada,tratamiento) |> 
  reframe(m = mean(peso < 4)*100,
          l = mean(peso >= 4 & peso < 6)*100,
          xl = mean(peso >= 6)*100,
          mediana = median(peso,na.rm=T)) |> 
  group_by(sitio,temporada) |> 
  mutate(dif_T0 = round(mediana-first(mediana),2))

diametro <- data_apariencia |> 
  group_by(sitio,temporada,tratamiento) |> 
  reframe(m = mean(diametro < 22)*100,
          l = mean(diametro >= 22 & peso < 24)*100,
          xl = mean(diametro >= 24)*100,
          mediana = median(diametro,na.rm=T)) |> 
  group_by(sitio,temporada) |> 
  mutate(dif_T0 = round(mediana-first(mediana),2))

color <- data_apariencia |> 
  group_by(sitio,temporada,tratamiento) |> 
  reframe(caoba = mean(between(color,4,5))*100,
          mediana = median(color,na.rm=T)) |> 
  group_by(sitio,temporada) |> 
  mutate(dif_T0 = round(mediana-first(mediana),2))

brix <- data_brix |> 
  group_by(sitio,temporada,tratamiento) |> 
  reframe(brix_optimo = mean(between(brix,19,20))*100,
          mediana = median(brix,na.rm=T)) |> 
  group_by(sitio,temporada) |> 
  mutate(dif_T0 = round(mediana-first(mediana),2))

daño <- data_daño |> 
  group_by(sitio,temporada,tratamiento) |> 
  reframe(mediana = median(daño,na.rm=T)) |> 
  group_by(sitio,temporada) |> 
  mutate(dif_T0 = round(mediana-first(mediana),2))

peso

diametro

color

brix

