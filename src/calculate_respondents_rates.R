setwd("/Users/hhakimjavadi/Dropbox (UFL)/FMT_microbiome_smichail/")
metadata<-read_excel("data/metadata/FMT Metadata 08.17.24.xlsx")

pucai= metadata %>% select (`Subject ID`, Group,`Treatment timepoint`, `Open label follow-up`, pucai, responder)

## pucai change in treatment group
pucai_tr=pucai %>% filter(Group == "Treatment")
pucai_tr_w=pucai_tr %>% pivot_wider( names_from = `Treatment timepoint`,
                                                                values_from = pucai,
                                                                values_fn = mean)
pucai_tr_w=pucai_tr_w %>% mutate(puci_change=Treatment_0wk-Treatment_4wk)

pucai_tr_w$new_respnd <- ifelse(
  (pucai_tr_w$puci_change >=20 | pucai_tr_w$Treatment_4wk<=10),
  "Yes",
  "No"
)

pucai_tr_w$new_respnd %>% table()




## pucai change in placebo group
pucai_pl=pucai %>% filter(Group == "Placebo")
pucai_pl_w=pucai_pl %>% pivot_wider( names_from = `Treatment timepoint`,
                                     values_from = pucai,
                                     values_fn = mean)
pucai_pl_w=pucai_pl_w %>% mutate(puci_change=Placebo_0wk-Placebo_4wk)

pucai_pl_w$new_respnd <- ifelse(
  (pucai_pl_w$puci_change >=20 | pucai_pl_w$Placebo_4wk<=10),
  "Yes",
  "No"
)

pucai_pl_w$new_respnd %>% table()



## pucai change in open label placebo-treatment group
pucai_pl_tr=pucai %>% filter(`Open label follow-up` == "Placebo-Treatment")
pucai_pl_tr_w=pucai_pl_tr %>% pivot_wider( names_from = `Treatment timepoint`,
                                     values_from = pucai,
                                     values_fn = mean)
pucai_pl_w=pucai_pl_w %>% mutate(puci_change=Placebo_0wk-Placebo_4wk)

pucai_pl_w$new_respnd <- ifelse(
  (pucai_pl_w$puci_change >=2-0 | pucai_pl_w$Placebo_4wk<=10),
  "Yes",
  "No"
)

pucai_pl_w$new_respnd %>% table()
