-- Tags common to both H.sapiens and S.scrofa genomes
create temp table c_tss AS (
  select tmp_hg95_ss9.* from tmp_hg95_ss9 inner join tmp_hg95_ncbi37
    on tmp_hg95_ss9.read_name = tmp_hg95_ncbi37.read_name
  );
--
select count_pattern(list_mism, '>') AS no_mism, 
  count(*) AS count_mism, 
  count(*)::numeric / (select count(*) from tmp_hg95_ss9) AS prop_mism
from tmp_hg95_ss9
group by no_mism
order by no_mism;

select count_pattern(list_mism, '>') AS no_mism, 
  count(*) AS count_mism, 
  count(*)::numeric / (select count(*) from tmp_hg95_ncbi37) AS prop_mism
from tmp_hg95_ncbi37
group by no_mism
order by no_mism;

select count_pattern(list_mism, '>') AS no_mism, 
  count(*) AS count_mism, 
  count(*)::numeric / (select count(*) from c_tss) AS prop_mism
from c_tss
group by no_mism
order by no_mism;
