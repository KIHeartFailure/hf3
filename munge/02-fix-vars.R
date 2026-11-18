escdata <- escdata %>%
  mutate(
    BASEDT = dmy(BASEDT),
    byear = factor(year(BASEDT)),
    DMVTYPE = fct_rev(fct_recode(DMVTYPE, AHF = "Acute HF event requiring iv therapy or mechanical support", CHF = "Out-patient visit")),
    DMPUA = fct_recode(DMPUA, "Significant primary valve disease" = "Significant Valve disease"),
    DCICCUTIME = if_else(DCICCU == 0, NA_real_, DCICCU),
    DCICCUYN = ynfac(if_else(DCICCU == 0, 0, 1)),
    DMHFHIST_DUR = factor(case_when(
      DMHFHIST == "No (now de Novo)" ~ "De novo",
      TRUE ~ str_remove(DMHFDUR, " months")
    ), levels = c("De novo", "< 6", "6-12", "> 12")),
    DMHFHOSP = ynfac(case_when(
      is.na(DMHFHIST) ~ NA_real_,
      DMHFHIST == "Yes, with previous hospitalization" ~ 1,
      TRUE ~ 0
    )),
    bsa1 = sqrt(DMHEIGHT * DMWEIGHT / 3600),
    bsa2 = sqrt(DMHEIGHT * DCWEIGHT / 3600),
    OPLAVOLI = OPLAVOL / bsa1,
    ITLAVOLI = ITLAVOL / bsa2,
    OPLVMI = OPLVM / bsa1,
    ITLVMI = ITLVM / bsa2,
    F1DIAB = factor(if_else(F1DIAB == "3", NA_character_, F1DIAB)),

    # Medications
    MDRASIARNIH = ynfac(case_when(
      is.na(MDACEH_C) | is.na(MDATH_C) | is.na(MDARNIH_C) ~ NA_real_,
      MDACEH_C == "Yes" | MDATH_C == "Yes" | MDARNIH_C == "Yes" ~ 1,
      TRUE ~ 0
    )),
    MDRASIARNID = ynfac(case_when(
      is.na(MDACED_C) | is.na(MDATD_C) | is.na(MDARNID_C) ~ NA_real_,
      MDACED_C == "Yes" | MDATD_C == "Yes" | MDARNID_C == "Yes" ~ 1,
      TRUE ~ 0
    )),
    MDSGLT2H = ynfac(case_when(
      is.na(MDNCV) ~ NA_real_,
      MDNCV == "No" ~ 0,
      is.na(MDDIORH) ~ NA_real_,
      MDDIORH == "SGLT2 (or 2/1) inhibitor" ~ 1,
      !is.na(MDDOSGH) ~ 1,
      TRUE ~ 0
    )),
    MDSGLT2D = ynfac(case_when(
      is.na(MDNCV) ~ NA_real_,
      MDNCV == "No" ~ 0,
      is.na(MDDIORD) ~ NA_real_,
      MDDIORD == "SGLT2 (or 2/1) inhibitor" ~ 1,
      !is.na(MDDOSGD) ~ 1,
      TRUE ~ 0
    )),
    MDGLP1H = ynfac(case_when(
      is.na(MDNCV) ~ NA_real_,
      MDNCV == "No" ~ 0,
      is.na(MDDNOGLH) ~ NA_real_,
      MDDNOGLH == "None" ~ 0,
      TRUE ~ 1
    )),
    MDGLP1D = ynfac(case_when(
      is.na(MDNCV) ~ NA_real_,
      MDNCV == "No" ~ 0,
      is.na(MDDNOGLD) ~ NA_real_,
      MDDNOGLD == "None" ~ 0,
      TRUE ~ 1
    )),
    all_vaso_na = if_all(all_of(c(
      "ITVASNO",
      "ITVASDOB",
      "ITVASDOP",
      "ITVASMIL",
      "ITVASENO",
      "ITVASLEV",
      "ITVASNOR",
      "ITVASEPI",
      "ITVASVAS",
      "ITVASOTH"
    )), is.na),

    # Outcomes
    startdtm = coalesce(DCDISCHDTD1, BASEDT, DMENROLDTD1),
    enddtm = pmin(F1DEATHDTD1, F1CONTDTD1, na.rm = T),
    time_death_fu = as.numeric(enddtm - startdtm),
    time_death_fu = if_else(time_death_fu <= 0 & DCVITAL == "Alive", 1, time_death_fu), 
    death_inhosp = ynfac(case_when(
      is.na(DMVTYPE) | DMVTYPE == "CHF" ~ NA_real_,
      DCVITAL == "Dead" ~ 1,
      DCVITAL == "Alive" | !is.na(F1VITAL) ~ 0
    )),
    death_fu = ynfac(case_when(
      F1DEATHDTD1 > F1CONTDTD1 ~ 0,
      F1VITAL == "Dead" ~ 1,
      F1VITAL == "Alive" ~ 0
    )),
    cvd_fu = case_when(
      is.na(death_fu) ~ NA_real_,
      death_fu == "No" | F1DTHCSP %in% c("Non-cardiovascular cause", "Systemic haemorrhage") ~ 0,
      TRUE ~ 1 # Unknown cause of death considered to be CVD
    ),
    cvdhfh_fu = case_when(
      is.na(F1HOSP) | is.na(cvd_fu) ~ NA_real_,
      TRUE ~ F1HOSP + cvd_fu
    ),
    cvd_fu = ynfac(cvd_fu)
  )

vasofunc <- function(data, vasovar) {
  data <- data %>%
    mutate(
      !!sym(vasovar) := ynfac(case_when(
        all_vaso_na ~ NA_real_,
        !!sym(vasovar) == "Yes" ~ 1,
        TRUE ~ 0
      ))
    )
}

escdata <- vasofunc(escdata, vasovar = "ITVASNO")
escdata <- vasofunc(escdata, vasovar = "ITVASDOB")
escdata <- vasofunc(escdata, vasovar = "ITVASDOP")
escdata <- vasofunc(escdata, vasovar = "ITVASMIL")
escdata <- vasofunc(escdata, vasovar = "ITVASENO")
escdata <- vasofunc(escdata, vasovar = "ITVASLEV")
escdata <- vasofunc(escdata, vasovar = "ITVASNOR")
escdata <- vasofunc(escdata, vasovar = "ITVASEPI")
escdata <- vasofunc(escdata, vasovar = "ITVASVAS")
escdata <- vasofunc(escdata, vasovar = "ITVASOTH")

diureticfunc <- function(data, diur1, diur2) {
  data <- data %>%
    mutate(
      tmp_furo1 = case_when(
        !!sym(diur1) == "Furosemide" ~ !!sym(paste0(diur1, "DO")),
        !!sym(diur1) == "Torasemide" ~ !!sym(paste0(diur1, "DO")) * 2,
        !!sym(diur1) == "Bumetanide" ~ !!sym(paste0(diur1, "DO")) * 40
      ),
      tmp_furo2 = case_when(
        !!sym(diur2) == "Furosemide" ~ !!sym(paste0(diur2, "DO")),
        !!sym(diur2) == "Torasemide" ~ !!sym(paste0(diur2, "DO")) * 2,
        !!sym(diur2) == "Bumetanide" ~ !!sym(paste0(diur2, "DO")) * 40
      ),
      !!sym(paste0(diur1, "DO_furo")) := if_else(
        rowSums(is.na(across(tmp_furo1:tmp_furo2))) == 2,
        NA_real_,
        rowSums(across(tmp_furo1:tmp_furo2), na.rm = TRUE)
      )
    ) %>%
    select(-tmp_furo1, -tmp_furo2)
  return(data)
}

escdata <- diureticfunc(escdata, diur1 = "MDDIURH", diur2 = "MDDIUR2H")
escdata <- diureticfunc(escdata, diur1 = "MDDIURD", diur2 = "MDDIUR2D")
escdata <- diureticfunc(escdata, diur1 = "F1DIURH", diur2 = "F1DIUR2H")
escdata <- diureticfunc(escdata, diur1 = "F1DIURD", diur2 = "F1DIUR2D")

labconvfunc <- function(var, unitconv, cf, add = 0) {
  escdata <- escdata %>%
    mutate(!!sym(var) := if_else(!!sym(paste0(var, "U")) == unitconv, !!sym(var) * cf + add, !!sym(var)))
}


escdata <- labconvfunc("OPWBC", "Cells/microL", 0.001)
escdata <- labconvfunc("OPCOL", "mg/dL", 0.02586)
escdata <- labconvfunc("OPHB", "g/dL", 10)
escdata <- labconvfunc("OPHB", "mmol/L", 16.114)
escdata <- labconvfunc("OPCOLLD", "mg/dL", 0.02586)
escdata <- labconvfunc("OPSERI", "µg/dL", 0.1791)
escdata <- labconvfunc("OPGL", "mg/dL", 1 / 18)
escdata <- labconvfunc("OPHB1", "mmol/mol", 0.0915, add = 2.15)
escdata <- labconvfunc("OPBNP", "mmol/mol", 3.4641)
escdata <- labconvfunc("OPCREA", "mg/dL", 88.42)
escdata <- labconvfunc("F1CREA", "mg/dL", 88.42)
escdata <- labconvfunc("OPNT", "mmol/mol", 8.457)

# BUN/UREA
escdata <- labconvfunc("OPURE", "mg/dL", 0.1665)
escdata <- labconvfunc("OPBUN", "mg/dL", 0.357)

escdata <- labconvfunc("OPURC", "mg/dL", 59.48)
escdata <- labconvfunc("OPTBIL", "mg/dL", 17.1)

escdata <- escdata %>%
  mutate(
    OPBUN = coalesce(OPBUN, OPURE),
    male = if_else(DMSEX == "Female", 0, 1),
    ethno = if_else(DMETHO == "Black", 1, 0),
    OPEGFR = nephro::CKDEpi.creat(OPCREA / 88.42, male, DMAGE, ethno),
    F1EGFR = nephro::CKDEpi.creat(F1CREA / 88.42, male, DMAGE, ethno)
  )

escdata <- escdata %>%
  mutate(across(where(is_character), factor)) %>%
  mutate(across(where(is.factor), droplevels))
