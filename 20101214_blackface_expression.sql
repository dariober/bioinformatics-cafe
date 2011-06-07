/*
Data for Figure 3 in Mol Ecol paper
*/

-- Log-fold change:

select distinct "LibraryID" from out_solexa_cumcountall_tpm;

-- Add group type to each library (C, R, S)
create temp table tag_expr AS(
    select out_solexa_cumcountall_tpm.*, "qrySolexaLambs"."TreatmentGroup"
    from out_solexa_cumcountall_tpm inner join "qrySolexaLambs" on
        out_solexa_cumcountall_tpm."LibraryID" = "qrySolexaLambs"."LibraryID"
    );
select * from tag_expr limit 10;

-- Average expression of each tag in each group:
drop table out_tag_treatment;
create table out_tag_treatment AS(
    select "TagChrPos", "TreatmentGroup", avg("CumCount"::double precision) AS avg_count
    from tag_expr
    group by "TreatmentGroup", "TagChrPos"
    -- Add RS group:
    UNION SELECT "TagChrPos", 'RS', avg("CumCount") AS avg_count
    from tag_expr
    where "TreatmentGroup" in ('R', 'S')
    group by "TagChrPos"
    );
select * from out_tag_treatment limit 1000;
select * from out_tag_treatment where "TagChrPos" like '10_16580536_1';

-- Cross tab: Tags as rows, TreatmentGroup as column, value is avg expr:
select cross_tab($$ select "TagChrPos", "TreatmentGroup", avg_count from out_tag_treatment $$, 'out_crosstab_tag_treatment');

select * from out_crosstab_tag_treatment limit 120;
select * from out_crosstab_tag_treatment where "TagChrPos" like '10_83553379_1';

-- Calculate log-fold changes and addE-values: (This is quite a badly designed query)
drop table logfc;
create temp table logfc AS(
  select 'CvsR' AS "TreatmentGroup", "SAGEoutput"."TagChrPos", log(2, "C"/"R") AS log2fc, "SAGEoutput"."E"
  from out_crosstab_tag_treatment inner join "SAGEoutput" on
       out_crosstab_tag_treatment."TagChrPos" = "SAGEoutput"."TagChrPos"
  where "SAGEoutput"."Dataset" = 'CvsR'

  union select 'CvsS' AS "TreatmentGroup", "SAGEoutput"."TagChrPos", log(2, "C"/"S") AS log2fc, "SAGEoutput"."E"
  from out_crosstab_tag_treatment inner join "SAGEoutput" on
       out_crosstab_tag_treatment."TagChrPos" = "SAGEoutput"."TagChrPos"
  where "SAGEoutput"."Dataset" = 'CvsS'

  union select 'RvsS' AS "TreatmentGroup", "SAGEoutput"."TagChrPos", log(2, "R"/"S") AS log2fc, "SAGEoutput"."E"
  from out_crosstab_tag_treatment inner join "SAGEoutput" on
       out_crosstab_tag_treatment."TagChrPos" = "SAGEoutput"."TagChrPos"
  where "SAGEoutput"."Dataset" = 'RvsS'
  );
copy logfc to 'D:/Tritume/logfc_evalue.csv' with csv header delimiter E'\t'; -- This has been moved to M:\Documents\Manuscripts\Molecular Ecology MS\Figure 3

select * from logfc limit 100;
select * from logfc where "TagChrPos" = '11_74134717_0'