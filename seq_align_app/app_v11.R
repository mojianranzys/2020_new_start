library(shiny)
library(stringr)
library(dplyr)
# UI
if(T){
  ui <- fluidPage(
    # App title ----
    titlePanel("Sequence Alignment "),
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        # Input: Select a file ----
        fileInput("file1", "Choose TXT File",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        actionButton("run", "Go to run"),
        # Horizontal line ----
        tags$hr(),
        
        # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),
        
        # Input: Select separator ----
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        # Input: Select quotes ----
        radioButtons("quote", "Quote",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = ""),
        
        numericInput(inputId = "obs", 
                     min = 3,
                     label = "Number of aa to view:", 
                     value = 4),
        downloadButton("downloadData", "Download ALL Align"),
        downloadButton("downloadData2", "Download Part Align"),
      ),
      mainPanel(
        h1("Amino Acid Sequence Alignment ", align="center",
           style ="font-family: 'times'; font-size: 22pt "),
        tabsetPanel(type = "tabs",
                    tabPanel("Summary",pre(tableOutput("summary"))),
                    tabPanel("sequence",pre(tableOutput("seq"))),
                    tabPanel("ALL Sequence Alignment Result",pre(tableOutput("all"))),
                    tabPanel("Strong The ALL Result",pre(tableOutput("Strong_all"))),
                    tabPanel("Part Sequence Alignment Result",pre(tableOutput("part"))),
                    tabPanel("Strong The Part Result",pre(tableOutput("Strong_part")))
        )
      )
    )
  )
}

# Server
server <- function(input, output) {
  d1 = eventReactive(input$run,{
    req(input$file1)
    AA_df <- read.table(input$file1$datapath,
                        header = input$header,
                        sep = input$sep,
                        #quote = input$quote,
                        stringsAsFactors=F,
                        check.names=F
    )
    return(AA_df)
  })
  #Step1:----------------- Sequence All Alignment ----------------
  # output1: view sequence file
  output$summary <- renderTable({
    dat1 = d1()
  })
  Seq_infor <- function(x){
    seq_valid <- NULL
    for(i in seq(nrow(x))){
      name <- paste0(x[i,1],"_",x[i,2])
      seq <- substr(x[i,4],5,75)
      seq_name <- cbind(name,seq)
      seq_valid <- rbind(seq_valid,seq_name)
    }
    rownames(seq_valid) <- seq_valid[,1]
    seq_valid <- seq_valid[,-1]
    seq_valid <- as.data.frame(seq_valid)
    return(seq_valid)
  }
  output$seq <- renderTable({
    dat1 = d1()
    dat2 = Seq_infor(dat1)
  },colnames = T,rownames = T)
  
  All_result <- function(x){
    #print(x)
    result <- NULL
    n <- input$obs-1
    if(nrow(x)>0){
      all_c <- NULL
      all_a <- NULL
      for(j in seq(str_length(x[1,])-n)){
        all_b <- substring(x[1,],j,j+n)
        all_a <- c(all_a,all_b)
      }
      all_a <- unique(all_a)
      #print(all_a)
      same_dat <- NULL
      for(i in seq(nrow(x)-1)){
        for(j in seq(length(all_a))){
          same_aa <- str_extract(x[i+1,],all_a[j])
          if(!is.na(same_aa)){
            same_dat <- c(same_dat,same_aa)
          }
        }
      }
      #print(length(same_dat))
      if(length(same_dat) == 0){
        all_align_result <- paste0("This Multiple peptide sequences All alignment dont have ",n+1," AA;")
      }
      if(length(same_dat) > 0){
        table_c <- as.data.frame(table(unlist(same_dat)))
        #print(table_c)
        num <- nrow(x)-1
        all_result <- table_c[table_c[,2]==num,]
        if(nrow(all_result) == 0){
          all_align_result <- paste0("This Multiple peptide sequences All alignment dont have ",n+1," AA;")
        }
        if(nrow(all_result)>0){
          for(i in seq(nrow(all_result))){
            result <- paste0(result,as.character(all_result[i,1]),sep=",")
          }
        }
        if(length(result) ==1){
          all_align_result <- paste0("The ",n+1," AA of this Multiple peptide sequences ALL alignment is:",result,";")
        }
      }
    }
    return(all_align_result)
  }
  
  # output: Sequence All alignment result:
  output$all <- renderTable({
    dat1 = d1()
    n = isolate(input$obs)-1
    dat2 = Seq_infor(dat1) 
    dat3 = All_result(dat2)
  })
  
  strong_all <- function(x){
    n <- input$obs-1
    if(nrow(x)>0){
      all_c <- NULL
      all_a <- NULL
      for(j in seq(str_length(x[1,])-n)){
        all_b <- substring(x[1,],j,j+n)
        all_a <- c(all_a,all_b)
      }
      all_a <- unique(all_a)
      #print(all_a)
      same_dat <- NULL
      for(i in seq(nrow(x)-1)){
        for(j in seq(length(all_a))){
          same_aa <- str_extract(x[i+1,],all_a[j])
          #print(same_aa)
          if(!is.na(same_aa)){
            same_dat <- c(same_dat,same_aa)
          }
        }
      }
      if(length(same_dat) == 0){
        name_all <- "NO"
        end_all_a <- "NO"
        all_b <- "NO"
        all_c <- "NO"
      }
      if(length(same_dat) > 0){
        table_c <- as.data.frame(table(unlist(same_dat)))
        num <- nrow(x)-1
        all_result <- table_c[table_c[,2]==num,]
        if(nrow(all_result) == 0){
          name_all <- "NO"
          end_all_a <- "NO"
          all_b <- "NO"
          all_c <- "NO"
        }
        if(nrow(all_result)>0){
          name_all <- NULL
          all_a <- NULL
          all_b <- NULL
          all_c <- NULL
          result <- NULL
          for(i in seq(nrow(all_result))){
            pos <- NULL
            for(j in seq(nrow(x))){
              aa_pos <- unlist(str_locate(x[j,],as.character(all_result[i,1])))
              pos <- rbind(pos,aa_pos)
            }
            #print(pos)
            dat <- NULL
            for(j in seq(nrow(pos))){
              if(!is.na(pos[j,1])){
                dat <- rbind(dat,as.character(x[j,]))
              }else{
                next
              }
            }
            pos_aa <- as.data.frame(na.omit(pos))
            print(pos_aa)
            max_start <- max(pos_aa[,1])
            min_end <- min(pos_aa[,2])
            #print(c(max_start,min_end))
            dat_t <- paste0("S_result>----------",as.character(all_result [i,1]),"----------")
            for(j in seq(nrow(dat))){
              if(pos_aa[j,1] <= max_start | pos_aa[j,2] >= min_end){
                len_start <- max_start - pos_aa[j,1]
                len_end <- pos_aa[j,2] - min_end
                add_start <- c(rep("#",len_start))
                add_end <- c(rep("#",len_end))
                dat_add <- c(add_start,dat[j,],add_end)
                dat_strc <- str_c(dat_add,collapse = "")
                #print(dat_strc)
                name <- rownames(which((gsub("#","",dat_strc))==x,arr.ind = T))
                dat_name <- paste0(name,">",dat_strc)
                #print(dat_name)
                dat_t <- rbind(dat_t,dat_name)
              }
            }
          #print(dat_t)
            name_seq <- as.data.frame(strsplit(dat_t,"[>]"))[1,]
            name_seq <- t(name_seq)
            rownames(name_seq) <- NULL
            name_all <- rbind(name_all,name_seq)
            #print(name_all)
            seq <- as.data.frame(strsplit(dat_t,"[>]"))[2,]
            seq <- t(seq)
            rownames(seq) <- NULL
            poly_aa <- unlist(strsplit(seq[1],"----------"))[2]
            dat_a <- as.data.frame(strsplit(sub(poly_aa," ",seq)," "))[1,]
            dat_a <- t(dat_a)
            rownames(dat_a) <- NULL
            b_rep <- length(seq)
            dat_b <- rep(poly_aa,b_rep)
            colnames(dat_b) <- NULL
            dat_c <- as.data.frame(strsplit(sub(poly_aa," ",seq)," "))[2,]
            dat_c <- t(dat_c)
            rownames(dat_c) <- NULL
            all_a <- rbind(all_a,dat_a)
            all_b <- c(all_b,dat_b)
            all_c <- rbind(all_c,dat_c)
          }
          end_all_a <- NULL
          max_a <- max(str_length(all_a[,1]))
          for(i in seq(length(all_a))){
            if((str_length(all_a[i])) <= max_a){ 
              add_a_len <- max_a - str_length(all_a[i])
              add_a <- rep("#",add_a_len)
              a <- str_c(c(add_a,all_a[i]),collapse = "")
              end_all_a <- rbind(end_all_a,a)
            }
          }
        }
      }
    }
    all_result <- cbind(name_all,end_all_a,all_b,all_c)
    all_result <- gsub("#"," ",all_result)
    return(all_result)
  }
  # output3:Sequence All Aligment result to strong
  output$Strong_all <- renderTable({
    dat1 = d1()
    n = isolate(input$obs)-1
    dat2 = Seq_infor(dat1) 
    dat3 = strong_all(dat2)
  },colnames = F)
  
  # download the sequence all align result
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$obs,"_polypeptide_all_alignment_result.txt")
    },
    content = function(file) {
      dat1 = d1()
      n = isolate(input$obs)-1
      dat2 = Seq_infor(dat1) 
      dat3 = strong_all(dat2)
      write.table(dat3,file,col.names = FALSE,row.names = FALSE)
    }
  )
  #Step2:----------------- Sequence Part Alignment ----------------
  Part_result <- function(x){
    result <- NULL
    n <- input$obs-1
    if(nrow(x)>0){
      part_c <- NULL
      for(i in seq(nrow(x))){
        part_a <- NULL
        for(j in seq(str_length(x[i,])-n)){
          part_b <- substring(x[i,],j,j+n)
          part_a <- c(part_a,part_b)
        }
        part_a <- unique(part_a)
        #print(part_a)
        same_dat <- NULL
        if(i <= (nrow(x)-1)){
          for(j in (i:(nrow(x)-1))){
            for(z in seq(length(part_a))){
              same_aa <- str_extract(x[j+1,],part_a[z])
              if(!is.na(same_aa)){
                same_dat <- c(same_dat,same_aa)
              }
            }
          }
        }else{
          next
        }
        part_c <- c(part_c,same_dat)
      }
      part_c <- unique(part_c)
      if(length(part_c) == 0){
        part_align_result <- paste0("This Multiple peptide sequences part alignment dont have ",n+1," AA;")
      }
      if(length(part_c)>0){
        for(i in seq(length(part_c))){
          result <- paste0(result,as.character(part_c[i]),sep=",")
        }
      }
      if(length(result) ==1){
        part_align_result <- paste0("The ",n+1," AA of this Multiple peptide sequences part alignment is:",result,";")
      }
      return(part_align_result)
    }
  }
  # output1:Sequence Part alignmnet result view:
  output$part <- renderTable({
    dat1 = d1()
    n = isolate(input$obs)-1
    dat2 = Seq_infor(dat1) 
    dat3 = Part_result(dat2)
  })
  
  strong_part <- function(x){
    n <- input$obs-1
    if(nrow(x)>0){
      part_c <- NULL
      for(i in seq(nrow(x))){
        part_a <- NULL
        for(j in seq(str_length(x[i,])-n)){
          part_b <- substring(x[i,],j,j+n)
          part_a <- c(part_a,part_b)
        }
        part_a <- unique(part_a)
        same_dat <- NULL
        if(i <= (nrow(x)-1)){
          for(j in (i:(nrow(x)-1))){
            for(z in seq(length(part_a))){
              same_aa <- str_extract(x[j+1,],part_a[z])
              if(!is.na(same_aa)){
                same_dat <- c(same_dat,same_aa)
              }
            }
          }
        }else{
          next
        }
        part_c <- c(part_c,same_dat)
      }
      #print(part_c)
      if(length(part_c) == 0){
        name_all <- "NO"
        end_all_a <- "NO"
        all_b <- "NO"
        all_c <- "NO"
      }
      if(length(part_c)>0){
        name_all <- NULL
        all_a <- NULL
        all_b <- NULL
        all_c <- NULL
        partc_table <- as.data.frame(table(part_c))
        result <- NULL
        for(i in seq(nrow(partc_table))){
          pos <- NULL
          for(j in seq(nrow(x))){
            aa_pos <- unlist(str_locate(x[j,],as.character(partc_table[i,1])))
            pos <- rbind(pos,aa_pos)
          }
          #print(pos)
          dat <- NULL
          for(j in seq(nrow(pos))){
            if(!is.na(pos[j,1])){
              dat <- rbind(dat,as.character(x[j,]))
            }else{
              next
            }
          }
          print(dat)
          pos_aa <- as.data.frame(na.omit(pos))
          max_start <- max(pos_aa[,1])
          min_end <- min(pos_aa[,2])
          dat_t <- paste0("S_result>----------",as.character(partc_table[i,1]),"----------")
          for(j in seq(nrow(dat))){
            if(pos_aa[j,1] <= max_start | pos_aa[j,2] >= min_end){
              len_start <- max_start - pos_aa[j,1]
              len_end <- pos_aa[j,2] - min_end
              add_start <- c(rep("#",len_start))
              add_end <- c(rep("#",len_end))
              dat_add <- c(add_start,dat[j,],add_end)
              dat_strc <- str_c(dat_add,collapse = "")
              #print(dat_strc)
              name <- rownames(which((gsub("#","",dat_strc))==x,arr.ind = T))
              dat_name <- paste0(name,">",dat_strc)
             # print(dat_name)
              dat_t <- rbind(dat_t,dat_name)
            }
          }
          #print(dat_t)
          name_seq <- as.data.frame(strsplit(dat_t,"[>]"))[1,]
          name_seq <- t(name_seq)
          rownames(name_seq) <- NULL
          name_all <- rbind(name_all,name_seq)
          seq <- as.data.frame(strsplit(dat_t,"[>]"))[2,]
          seq <- t(seq)
          rownames(seq) <- NULL
          poly_aa <- unlist(strsplit(seq[1],"----------"))[2]
          dat_a <- as.data.frame(strsplit(sub(poly_aa," ",seq)," "))[1,]
          dat_a <- t(dat_a)
          rownames(dat_a) <- NULL
          b_rep <- length(seq)
          dat_b <- rep(poly_aa,b_rep)
          colnames(dat_b) <- NULL
          dat_c <- as.data.frame(strsplit(sub(poly_aa," ",seq)," "))[2,]
          dat_c <- t(dat_c)
          rownames(dat_c) <- NULL
          all_a <- rbind(all_a,dat_a)
          all_b <- c(all_b,dat_b)
          all_c <- rbind(all_c,dat_c)
        }
        end_all_a <- NULL
        max_a <- max(str_length(all_a[,1]))
        for(i in seq(length(all_a))){
          if((str_length(all_a[i])) <= max_a){ 
            add_a_len <- max_a - str_length(all_a[i])
            add_a <- rep("#",add_a_len)
            a <- str_c(c(add_a,all_a[i]),collapse = "")
            end_all_a <- rbind(end_all_a,a)
          }
        }
      }
    }
    all_result <- cbind(name_all,end_all_a,all_b,all_c)
    all_result <- gsub("#"," ",all_result)
    return(all_result)
  }
  #output2:Sequence Part alignmnet result to strong view:
  output$Strong_part <- renderTable({
    dat1 = d1()
    n = isolate(input$obs)-1
    dat2 = Seq_infor(dat1) 
    dat3 = strong_part(dat2)
  },colnames = F)
  
  # download the sequence Part align result
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0(input$obs,"_polypeptide_part_alignment_result.txt")
    },
    content = function(file) {
      dat1 = d1()
      n = isolate(input$obs)-1
      dat2 = Seq_infor(dat1)
      dat3 = strong_part(dat2)
      write.table(dat3,file,col.names = FALSE,row.names = FALSE)
    }
  )
}
shinyApp(ui,server)
