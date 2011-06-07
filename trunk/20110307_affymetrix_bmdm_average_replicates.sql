/*
  Average affymetrix arrays to have the mean of the technical replicates.
  NOTE: s024_LWL5_24h_repl_2 excluded.
*/
drop table if exists affyarray_bmdm_avg;
create table affymetrix.affyarray_bmdm_avg AS(
select 
    pig || '_' || time_point || 'h'::text AS array_id,
    pig::text, 
    time_point::int,
    probe_set_id,
    avg( log_intensity ) AS log2_intensity
from 
    ( select *, log(2.0, intensity::numeric)::numeric AS log_intensity from affyarray_bmdm ) as t1
inner join affyarray_design on affyarray_design.slide_id = t1.slide_id
where affyarray_design.slide_id != 's024_LWL5_24h_repl_2'
group by pig, time_point, probe_set_id
order by pig, time_point, probe_set_id
);

comment on table affyarray_bmdm_avg is 'Expression data for the bmdm affymetrix experiment with technical replicates log2 transformed and then averaged. Original data is in table affyarray_bmdm. Array s024_LWL5_24h_repl_2 excluded. See 20110307_affymetrix_bmb_average_replicates.sql'
alter table affyarray_bmdm_avg add constraint pk_probes primary key (probe_set_id, time_point, pig);
vacuum analyze affyarray_bmdm_avg;
