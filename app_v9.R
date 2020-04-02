library(shiny)
library(stringr)
library(dplyr)

# Define UI for data upload app ----
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
                     selected = "\t"),
        
        # Input: Select quotes ----
        radioButtons("quote", "Quote",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = "'"),
        
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
                    #tabPanel("Split to polypeptide", pre(tableOutput("split"))),
                    tabPanel("ALL Sequence Alignment Result",pre(tableOutput("all"))),
                    #tabPanel("Line of The ALL Result",pre(tableOutput("line_all"))),
                    tabPanel("Strong The ALL Result",pre(tableOutput("Strong_all"))),
                    tabPanel("Part Sequence Alignment Result",pre(tableOutput("part"))),
                    #tabPanel("Line of The Part Result",pre(tableOutput("line_part"))),
                    tabPanel("Strong The Part Result",pre(tableOutput("Strong_part")))
        )
      )
    )
  )
}

server <- function(input, output) {
  d1 = eventReactive(input$run,{
    req(input$file1)
    AA_df <- read.table(input$file1$datapath,
                        header = input$header,
                        sep = input$sep,
                        quote = input$quote,
                        stringsAsFactors=F
    )
    rownames(AA_df)=paste('Sequence_',rownames(AA_df),sep='') 
    return(AA_df)
  })
  # output1: ---
  output$summary <- renderTable({
    dat1 = d1()
  },rownames = T,colnames = F)
  
  # step0: split string_aa to char
  split_str <- function(x){
    len_aa <- NULL
    for(i in seq(nrow(x))){
      raw_split <- unlist(strsplit(x[i,],split = ""))
      len_aa <- rbind(len_aa,length(raw_split))
    }
    #print(max(len_aa))
    raw_dat <- NULL
    for(i in seq(nrow(x))){
      raw_split <- unlist(strsplit(x[i,],split = ""))
      #  print(raw_split)
      if(length(raw_split) <- max(len_aa)){
        rep_n <- max(len_aa) - length(raw_split)
        raw_split <- t(cbind(raw_split,rep("#",rep_n)))
        raw_dat <- rbind(raw_dat,raw_split)
      }
    }
    #print(raw_dat)
    colnames(raw_dat) <- NULL
    rownames(raw_dat) <- NULL
    raw_dat[which(is.na(raw_dat),arr.ind = T)] <- "#"
    #    print(length(raw_dat))
    return(raw_dat)
  }
  
  # step1:split a long string to 4~10 aa
  split_fun <- function(x){
    dat <- NULL
    for(j in seq(length(x))){
      n <- input$obs-1
      m <- j+n
      aa <- str_c(x[j:m],collapse = "")
      dat <- cbind(dat,aa)
    }
    colnames(dat) <- NULL
    dat <- as.data.frame(dat)
    return(dat)
  }
  poly <- function(x){
    polypeptide_data <- NULL
    for(i in seq(nrow(x))){
      polypeptide_data <- rbind(polypeptide_data,split_fun(x[i,]))
      #      print(polypeptide_data)
    }
    polypeptide_data <- t(na.omit(t(polypeptide_data)))
    return(polypeptide_data)
  }
  # output2: --- NO View
  if(F){
    output$split <- renderTable({
      dat1 = d1()
      dat2 = split_str(dat1)
      dat3 = split_fun(dat2)
      dat4 = poly(dat2)
    })
  }
  # step2 :find Same continuous polypeptide
  find_same <- function(x){
    n <- input$obs-1
    AA_1_df <- NULL
    for(i in seq(nchar(x[1,])-n)){
      AA_1 <- substring(x[1,],i,i+n)
      AA_1_df <- c(AA_1_df,AA_1)
    }
    AA_1_df <- unique(AA_1_df)
    same_dat <- NULL
    for(i in seq(nrow(x)-1)){
      for(j in seq(length(AA_1_df))){
        same_aa <- str_extract(x[i+1,],AA_1_df[j])
        if(!is.na(same_aa)){
          same_dat <- c(same_dat,same_aa)
        }
      }
    }
    return(same_dat)
  }
  
  #step3 polypeptide Sequence Alignment
  Sequence_Alignment <- function(x,y,n){
    polypeptide_data <- NULL
    AA_same <- find_same(y)
    #print(length(AA_same))
    if(length(AA_same)==0){
      Align_result <- paste0("This Multiple peptide sequences All alignment dont have ",n+1," AA;")
    }
    if(length(AA_same)>=1){
      same_aa_table <- as.data.frame(table(unlist(AA_same)))
      same_aa_table <- same_aa_table[!(grepl("#",same_aa_table[,1])),]
      result <- NULL
      num_data <- nrow(x)
      num <- num_data-1
      bb <- same_aa_table[same_aa_table$Freq == num,]
      if(nrow(bb) ==0){
        Align_result <- paste0("This Multiple peptide sequences All alignment dont have ",n+1," AA;")
      }
      if(nrow(bb)>0){
        for(i in seq(nrow(same_aa_table))){
          if(same_aa_table$Freq[i] >= num){
            result <- paste0(result,as.character(same_aa_table$Var1[i]),sep=",")
          }
        }
      }
      if(length(result) ==1){
        Align_result <- paste0("The ",n+1," AA of this Multiple peptide sequences All alignment is:",result,";")
      }
    }
    #    print(Align_result)
    return(Align_result)
  }
  output$all <- renderTable({
    dat1 = d1()
    dat2 = split_str(dat1)
    n = isolate(input$obs)-1
    Sequence_Alignment(dat2,dat1,n)
  })
  # view sequence All alignment result
  Seq_align_view <- function(x,y){
    dat_all <- NULL
    #  dat_t <- NULL
    judge <- grep(":",x)
    if(length(judge) >0){
      t <- unlist(strsplit(unlist(strsplit(gsub(";","",x),"[:]"))[2],","))
      #print(t)
      AA_all <- NULL
      for(i in seq(nrow(y))){
        df1 <- str_c(y[i,],collapse = "")
        AA_all <- rbind(AA_all,df1)
      }
      #print(AA_all)
      for(i in seq(length(t))){
        dat_t <- paste0("----------",t[i],"----------")
        pos <- NULL
        for(j in seq(nrow(AA_all))){
          pos_t <- unlist(str_locate(AA_all[j,],t[i]))
          pos <- rbind(pos,pos_t)
        }
        #print(pos)
        colnames(pos) <- c("start","end")
        max_start <- max(pos[,1])
        min_end <- min(pos[,2])
        for(j in seq(nrow(y))){
          if(pos[j,1] < max_start | pos[j,2] > min_end){
            len_start <- max_start - pos[j,1]
            len_end <- pos[j,2] - min_end
            add_start <- c(rep("#",len_start))
            add_end <- c(rep("#",len_end))
            dat_add <- c(add_start,y[j,],add_end)
            dat_strc <- str_c(dat_add,collapse = "")
            dat_t <- rbind(dat_t,dat_strc)
          }
        }
        dat_all <- rbind(dat_all,dat_t)
        dat_result <- gsub("#"," ",dat_all)
      }
      return(dat_result)
    }
    if(length(judge) == 0){
      dat_result <- paste("The All Sequence alingment is NOthing To View")
      return(dat_result)
    }
  }
  # output NO View
  if(F){
    output$line_all <- renderTable({
      dat1 = d1()
      dat2 = split_str(dat1)
      n = isolate(input$obs)-1
      dat3 = Sequence_Alignment(dat2,dat1,n)
      dat4 = Seq_align_view(dat3,dat2)
    },colnames = F,fileEncoding = "utf-8")
    
  }
  
  
  # split the all align result to part 3,this is 1:
  Color_result_a <- function(x,y){
    dat_a  <- NULL
    judge <- grep(":",x)
    if(length(judge) >0){
      t <- unlist(strsplit(unlist(strsplit(gsub(";","",x),"[:]"))[2],","))
      AA_all <- NULL
      for(i in seq(nrow(y))){
        df1 <- str_c(y[i,],collapse = "")
        AA_all <- rbind(AA_all,df1)
      }
      #print(AA_all)
      dat_all <- NULL
      for(i in seq(length(t))){
        dat_t <- paste0("----------",t[i],"----------")
        pos <- NULL
        for(j in seq(nrow(AA_all))){
          pos_t <- unlist(str_locate(AA_all[j,],t[i]))
          pos <- rbind(pos,pos_t)
        }
        colnames(pos) <- c("start","end")
        max_start <- max(pos[,1])
        min_end <- min(pos[,2])
        for(j in seq(nrow(y))){
          if(pos[j,1] < max_start | pos[j,2] > min_end){
            len_start <- max_start - pos[j,1]
            len_end <- pos[j,2] - min_end
            add_start <- c(rep("#",len_start))
            add_end <- c(rep("#",len_end))
            dat_add <- c(add_start,y[j,],add_end)
            dat_strc <- str_c(dat_add,collapse = "")
            dat_t <- rbind(dat_t,dat_strc)
          }
        }
        #print(dat_t)
        color_a <- as.data.frame(strsplit(sub(t[i]," ",dat_t)," "))[1,]
        color_a <- t(color_a)
        rownames(color_a) <- NULL
        dat_a <- rbind(dat_a,color_a)
      }        
      end_all_a <- NULL
      max_a <- max(str_length(dat_a[,1]))
      for(i in seq(length(dat_a))){
        if((str_length(dat_a[i])) <= max_a){
          add_a_len <- max_a - str_length(dat_a[i])
          add_a <- rep("#",add_a_len)
          a <- str_c(c(add_a,dat_a[i]),collapse = "")
          end_all_a <- rbind(end_all_a,a)
        }
      }
      end_all_a <- gsub("#"," ",end_all_a)
      #end_all_a <- gsub("Polypeptide"," ",end_all_a)
      return(end_all_a)
    }
    if(length(judge) == 0){
      dat_a <- paste("NO")
      return(dat_a)
    }
  }
  # This is 2:
  Color_result_b <- function(x,y){
    m <- nrow(y)+1
    #print(m)
    judge <- grep(":",x)
    dat_b <- NULL
    if(length(judge) >0){
      t <- unlist(strsplit(unlist(strsplit(gsub(";","",x),"[:]"))[2],","))
      for(i in seq(length(t))){
        color_b <- rep(t[i],m)
        colnames(color_b) <- NULL
        #print(color_b)
        dat_b <- c(dat_b,color_b)
      }
      return(dat_b)
      
    }
    if(length(judge) == 0){
      dat_b <- paste("NO")
      return(dat_b)
    }
  }
  # This is 3:
  Color_result_c <- function(x,y){
    dat_c <- NULL
    judge <- grep(":",x)
    if(length(judge) >0){
      t <- unlist(strsplit(unlist(strsplit(gsub(";","",x),"[:]"))[2],","))
      AA_all <- NULL
      for(i in seq(nrow(y))){
        df1 <- str_c(y[i,],collapse = "")
        AA_all <- rbind(AA_all,df1)
      }
      #print(AA_all)
      for(i in seq(length(t))){
        dat_t <- paste0("----------",t[i],"----------")
        pos <- NULL
        for(j in seq(nrow(AA_all))){
          pos_t <- unlist(str_locate(AA_all[j,],t[i]))
          pos <- rbind(pos,pos_t)
        }
        colnames(pos) <- c("start","end")
        max_start <- max(pos[,1])
        min_end <- min(pos[,2])
        for(j in seq(nrow(y))){
          if(pos[j,1] < max_start | pos[j,2] > min_end){
            len_start <- max_start - pos[j,1]
            len_end <- pos[j,2] - min_end
            add_start <- c(rep("#",len_start))
            add_end <- c(rep("#",len_end))
            dat_add <- c(add_start,y[j,],add_end)
            dat_strc <- str_c(dat_add,collapse = "")
            dat_t <- rbind(dat_t,dat_strc)
          }
        }
        color_c <- as.data.frame(strsplit(sub(t[i]," ",dat_t)," "))[2,]
        color_c <- t(color_c)
        rownames(color_c) <- NULL
        dat_c <- rbind(dat_c,color_c)
      }
      dat_c <- gsub("#"," ",dat_c)
      return(dat_c)
    }
    if(length(judge) == 0){
      dat_c <- paste("NO")
      return(dat_c)
    }
  }
  output$Strong_all <- renderTable({
    dat1 = d1()
    dat2 = split_str(dat1)
    n = isolate(input$obs)-1
    dat3 = Sequence_Alignment(dat2,dat1,n)
    dat4 = Seq_align_view(dat3,dat2)
    dat5 = Color_result_a(dat3,dat2)
    dat6 = Color_result_b(dat3,dat2)
    dat7 = Color_result_c(dat3,dat2)
    dat8 = as.data.frame(cbind(dat5,dat6,dat7))
    #dat9 = gsub("#"," ",dat8)
  },colnames = F)
  # Two:Sequence part Alignmnet
  # step1: view the sequence part alignmnet result
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
  output$part <- renderTable({
    dat1 = d1()
    n = isolate(input$obs)-1
    dat2 = Part_result(dat1)
  })
  # Line the result:
  Part_sequence_align <- function(x){
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
      if(length(part_c) == 0){
        result <- paste("The Part sequence alignment result is NOthing To View ")
      }
      if(length(part_c) > 0){
        partc_table <- as.data.frame(table(part_c))
        for(i in seq(nrow(partc_table))){
          pos <- NULL
          for(j in seq(nrow(x))){
            aa_pos <- unlist(str_locate(x[j,],as.character(partc_table[i,1])))
            pos <- rbind(pos,aa_pos)
          }
          dat <- NULL
          for(j in seq(nrow(pos))){
            if(!is.na(pos[j,1])){
              dat <- rbind(dat,x[j,])
            }else{
              next
            }
          }
          pos_aa <- as.data.frame(na.omit(pos))
          #print(pos_aa)
          max_start <- max(pos_aa[,1])
          min_end <- min(pos_aa[,2])
          #print(c(max_start,min_end))
          result <- NULL
          dat_t <- paste0("----------",as.character(partc_table[i,1]),"----------")
          for(j in seq(nrow(dat))){
            if(pos_aa[j,1] < max_start | pos_aa[j,2] > min_end){
              len_start <- max_start - pos_aa[j,1]
              len_end <- pos_aa[j,2] - min_end
              add_start <- c(rep("#",len_start))
              add_end <- c(rep("#",len_end))
              dat_add <- c(add_start,dat[j,],add_end)
              dat_strc <- str_c(dat_add,collapse = "")
              dat_t <- rbind(dat_t,dat_strc)
              result <- rbind(result,dat_t)
            }
          }
        }
        result <- gsub("#"," ",dat_t)
      }
    }
    return(result)
  }
  # output NO view
  if(F){
    output$line_part <- renderTable({
      dat1 = d1()
      n = isolate(input$obs)-1
      dat2 = Part_sequence_align(dat1)
    },colnames = F,fileEncoding = "utf-8")
  }
  # split the sequence part align result to 3 part,and combine it:
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
              dat <- rbind(dat,x[j,])
            }else{
              next
            }
          }
          pos_aa <- as.data.frame(na.omit(pos))
          #print(pos_aa)
          max_start <- max(pos_aa[,1])
          min_end <- min(pos_aa[,2])
          #print(c(max_start,min_end))
          dat_t <- paste0("----------",as.character(partc_table[i,1]),"----------")
          for(j in seq(nrow(dat))){
            if(pos_aa[j,1] < max_start | pos_aa[j,2] > min_end){
              len_start <- max_start - pos_aa[j,1]
              len_end <- pos_aa[j,2] - min_end
              add_start <- c(rep("#",len_start))
              add_end <- c(rep("#",len_end))
              dat_add <- c(add_start,dat[j,],add_end)
              dat_strc <- str_c(dat_add,collapse = "")
              name <- rownames(which((gsub("#","",dat_strc))==x,arr.ind = T))
              dat_name <- paste0(name,">",dat_strc)
              dat_t <- rbind(dat_t,dat_name)
            }
          }
          print(dat_t)
          name_seq <- as.data.frame(strsplit(dat_t,"[>]"))[1,]
          name_seq <- t(name_seq)
          rownames(name_seq) <- NULL
          name_seq <- gsub(name_seq[grep("-",name_seq),],"",name_seq)
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
    all_result <- cbind(name_all,end_all_a,all_b,all_c)
    all_result <- gsub("#"," ",all_result)
    #all_result <- gsub("Polypeptide"," ",all_result)
    return(all_result)
  }
  
  output$Strong_part <- renderTable({
    dat1 = d1()
    n = isolate(input$obs)-1
    dat2 = strong_part(dat1)
  },colnames = F)
  # download the sequence all align result
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$obs,"_polypeptide_all_alignment_result.txt")
    },
    content = function(file) {
      dat1 = d1()
      dat2 = split_str(dat1)
      n = isolate(input$obs)-1
      dat3 = Sequence_Alignment(dat2,dat1,n)
      dat4 = Seq_align_view(dat3,dat2)
      dat5 = Color_result_a(dat3,dat2)
      dat6 = Color_result_b(dat3,dat2)
      dat7 = Color_result_c(dat3,dat2)
      dat8 = as.data.frame(cbind(dat5,dat6,dat7))
      write.table(dat8,file,col.names = FALSE,row.names = FALSE)
    }
  )
  # download the sequence Part align result
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0(input$obs,"_polypeptide_part_alignment_result.txt")
    },
    content = function(file) {
      dat1 = d1()
      n = isolate(input$obs)-1
      dat2 = strong_part(dat1)
      write.table(dat2,file,col.names = FALSE,row.names = FALSE)
    }
  )
}
shinyApp(ui,server)
