
###################################################
### Converrting a SAM result into a tibble form ###
###################################################

convert_vector <- function(vector,name){
  vector %>%
    as_tibble %>%
    mutate(year = as.integer(names(vector))) %>%
    mutate(type="VPA",sim="s0",stat=name,age=NA)
}

# VPAや将来予測の結果をtidyに ----
convert_df <- function(df,name){
  df %>%
    as_tibble %>%
    mutate(age = as.numeric(rownames(df))) %>%
    gather(key=year, value=value, -age, convert=TRUE) %>%
    group_by(year) %>%
    #        summarise(value=sum(value)) %>%
    mutate(type="VPA",sim="s0",stat=name)
}

#' SAMの結果オブジェクトをtibble形式に変換する関数
#'
#' @param vpares vpaの結果のオブジェクト
#' @encoding UTF-8
#'
#'
#' @export
#'
convert_sam_tibble <- function(samres) {

  if (class(samres) == "vpa") {
    total.catch <- colSums(samres$input$dat$caa*samres$input$dat$waa,na.rm=T)
  } else{
    total.catch <- colSums(samres$caa*samres$input$dat$waa,na.rm=T)
  }
  U <- total.catch/colSums(samres$baa, na.rm=T)

  SSB <- convert_vector(colSums(samres$ssb,na.rm=T),"SSB") %>%
    dplyr::filter(value>0&!is.na(value))
  Biomass <- convert_vector(colSums(samres$baa,na.rm=T),"biomass") %>%
    dplyr::filter(value>0&!is.na(value))
  FAA <- convert_df(samres$faa,"fishing_mortality") %>%
    dplyr::filter(value>0&!is.na(value))
  Recruitment <- convert_vector(colSums(samres$naa[1,,drop=F]),"Recruitment") %>%
    dplyr::filter(value>0&!is.na(value))

  Fratio <- NULL

  all_table <- bind_rows(SSB,
                         Biomass,
                         convert_vector(U[U>0],"U"),
                         convert_vector(total.catch[total.catch>0],"catch"),
                         convert_df(samres$naa,"fish_number"),
                         FAA,
                         convert_df(samres$input$dat$waa,"weight"),
                         convert_df(samres$input$dat$maa,"maturity"),
                         convert_df(samres$input$dat$caa,"catch_number"),
                         convert_df(samres$input$dat$M,  "natural_mortality"),
                         Recruitment,
                         Fratio) %>%
    mutate(type = "SAM")
}

#
# retro_plot = function(data,start_year=2005,scale=1000,mohn_res=NULL) {
#   data2 = data %>% dplyr::filter(stat %in% c("SSB","biomass","Recruitment","fishing_mortality")) %>%
#     group_by(id,year,stat) %>%
#     summarise(value=mean(value)) %>%
#     ungroup() %>%
#     mutate(value = if_else(stat == "fishing_mortality",value,value/1000)) %>%
#     mutate(terminal_year = max(data$year)-id) %>%
#     filter(year <= terminal_year) %>%
#     mutate(stat_f = factor(stat,levels=c("biomass","SSB","Recruitment","fishing_mortality"),labels=c("Biomass","SSB","Recruitment","F"))) %>%
#     filter(year >= start_year) %>%
#     filter(year <= as.numeric(terminal_year)) %>%
#     mutate(terminal_year = factor(terminal_year, levels=unique(terminal_year)))
#   g1 = ggplot(data=data2,aes(x=year,y=value))+
#     geom_path(aes(colour=terminal_year,linetype=terminal_year),size=1)+
#     facet_wrap(vars(stat_f),scales="free_y",ncol=2)+
#     theme_SH()+theme_bw(base_size=14)+theme(legend.position="none")+
#     xlab("Year") + ylab("")+
#     ylim(0,NA)
#   if (!is.null(mohn_res)) {
#     mohn = tibble(rho = mohn_res, stat = names(mohn_res)) %>%
#       dplyr::filter(stat != "N") %>%
#       mutate(stat_f = factor(stat,levels=c("B","SSB","R","F"),labels=c("Biomass","SSB","Recruitment","F"))) %>%
#       mutate(label=sprintf("rho == %.2f",rho)) %>%
#       full_join(data2 %>% group_by(stat_f) %>% summarise(ymax = max(value)))
#     g1 = g1 + geom_text(data=mohn,parse=TRUE,aes(x=start_year,y=ymax,label=label,hjust=0,vjust=1))
#   }
#   g1
# }
