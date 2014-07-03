package aligner;

import static org.junit.Assert.*;

import java.io.FileNotFoundException;

import org.junit.Test;


public class FKModTest {

	Aligner aln= new Aligner();
	String matrixFile= aln.getTmpFileFromResourceFile("nuc_N0.txt");
	FKMod fkmod= new FKMod();
	String query;
	String ref;
	
	@Test
	public void getBarcodesFromAlignment() throws FileNotFoundException{
		
		query= "GGGAAGCAAAATCGAGAATCCCGGTGCCGATACCCACTCTTGCAGAATCTATCACAACGTGTCAGATATATACATCCGAT";
		ref=   "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";
		aln.align(query, ref, matrixFile);
	
		String firstBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_FIRST_BARCODE);
		
		assertEquals("GGGAAGCAAA", firstBarcode);

		// With gaps in the query:
		query= "AAGCAAAATCGAGAATCCCGGTGCCGATACCCACTCTTGCAGAATCTATCACAACGTGTCAGATATATACATCCGAT";
		ref=   "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";
		aln.align(query, ref, matrixFile);
		firstBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_FIRST_BARCODE);
		assertEquals("---AAGCAAA", firstBarcode);
		
		// With gaps in the reference:
		query= "GGGAAGCAAATTTATCGAGAATCCCGGTGCCGATACCCACTCTTGCAGAATCTATCACAACGTGTCAGATATATACATCCGAT";
		ref=      "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";
		aln.align(query, ref, matrixFile);
		firstBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_FIRST_BARCODE);
		assertEquals("GGGAAGCAAA", firstBarcode); // Note that the gap could be at the beginning of the read! In fact it would be preferable

		// Second barcode truncated:
		// NB: Truncated barcodes should not be used as the aligner will probably split it.
		// E.g.:
		// Q: AAGCAAATTTATCGAGAATCCCGGTGCCGATACCCACTCTTGCAGAATCT-----------------------ATC----
		// R: NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT
		query=    "AAGCAAATTTATCGAGAATCCCGGTGCCGATACCCACTCTTGCAGAATCTATC";
		ref=      "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";
		aln.align(query, ref, matrixFile);
		String secondBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_SECOND_BARCODE);
		assertEquals("TCT-------", secondBarcode);
		
		// Revcomp
		query= "ATCGGATGTATATATCTGACACGCACAAACCTCTTCTGCAAGAGTGGGTATCGGCACCGGGATTCTCGATTCCGACCCC";
		ref=   "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";
		aln.align(query, ref, matrixFile);
		firstBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_FIRST_BARCODE);
		assertEquals("-GGGGTCGGA", firstBarcode);
	}
	
	@Test
	public void getModifiedBases() throws FileNotFoundException{
	
		query= "GGGAAGCAAAATCGAGAATCCCGGTGCCGATACCCACTCTTGTAGAATCTATCACAACGTGTCAGATATATACATCCGAT";
		ref=   "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";
		aln.align(query, ref, matrixFile);
		String modifiedBases= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_MODIFICATIONS);
		assertEquals("CT", modifiedBases);

		query= "GGGAAGCAAAATCGAGAATCCCGGTGCCGATACCCACT";
		ref=   "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";
		aln.align(query, ref, matrixFile);
		modifiedBases= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_MODIFICATIONS);
		assertEquals("C-", modifiedBases);

		// Rev comp
		// NB: Aligned to rev strand will return bases reverse complemented. So if the read is
		// TC you get GA
		query= "ATCGGATGTATATATCTGACACGNNNNNNNNNNTTCTtCAAGAGTcGGTATCGGCACCGGGATTCTCGATNNNNNNNNNN";
		ref=   "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";
		aln.align(query, ref, matrixFile);
		modifiedBases= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_MODIFICATIONS);
		assertEquals("ga", modifiedBases);
		
	}
	
	
}
