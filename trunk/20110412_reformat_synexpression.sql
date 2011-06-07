/*
  Convert table with probes on the the same row (as string) to one line per probe
*/

-- Table to convert
select * from attract_synexpression limit 100;

-- Convert string to array
drop table if exists syn;
create temp table syn AS(
    select 
       keggid,        -- row(s) to be repated along the array elements
       pway,
       pathway,
       pway_pvalue,
       cluster_group,
       string_to_array(syn_groups, ', ') AS probes
    from attract_synexpression
);
select * from syn limit 10;

-- Convert array to rows
-- drop table affymetrix.attract_synexpression_probes;
CREATE TABLE affymetrix.attract_synexpression_probes AS(
    select 
        keggid,        -- row(s) to be repated along the array elements
        pway,
        pathway,
        pway_pvalue,
        cluster_group,
        probes[i] AS probe_set_id -- Array whose elements will be rows (note [i])
    from syn, -- Input table
    -- generate_series() produce the series of indexes for the array elements.
    generate_series(1, 
         -- This query gets the length of the longest array in the table
        ( select max(array_upper(probes, 1)) from syn) ) i 
    -- Discard rows where the array index is larger than the array.
    where probes[i] is not null
    order by keggid, probe_set_id
);
create unique index indx_synexprs on attract_synexpression_probes (probe_set_id, cluster_group, keggid);
comment on table attract_synexpression_probes is 'Output of R package attract, function findSynexprs. Probe clusters within signficantly enriched KEGG pathways. See 20110413_attract_pig_affymetrix.R. This table is the same as attract_synexpression reformatted using 20110412_reformat_synexpression.sql'

select * from attract_synexpression_probes limit 10;
