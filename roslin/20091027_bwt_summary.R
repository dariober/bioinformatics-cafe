library(RODBC)

con<- odbcConnect("pgCAGE", case="nochange")

no_hits<- sqlQuery(con, "
  create temp table no_hits  AS (
    select job_id, read_name, count(*) AS no_hits
    from bowtie 
    group by job_id, read_name
  );
  " )

no_hits2<- sqlQuery(con, "  
  select job_id, no_hits, count(*) AS count_hits
  from no_hits
  group by job_id, no_hits;
  ")

hit_mat<- matrix(no_hits2$count_hits, ncol= 2)

tot_hit<- apply(hit_mat, 2, sum)

prop_hit_mat<- prop.table(hit_mat, 2)

par(mfrow= c(1,2))
barplot(t(hit_mat), beside= TRUE, names.arg= 1:nrow(hit_mat), 
  xlab= "Number of hits", ylab= "No. mapped tags", col= c("grey90", "grey30"))
legend("topright", legend= c("hg95ctrls vs S. scrofa", "hg95ctrls vs H. sapiens"), 
  fill= c("grey90", "grey30"))
barplot(t(prop_hit_mat), beside= TRUE, names.arg= 1:nrow(hit_mat), 
  xlab= "Number of hits", ylab= "Proportion of mapped tags", 
  col= c("grey90", "grey30"))  

odbcClose(con)

