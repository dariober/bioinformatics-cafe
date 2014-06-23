package bisReadBias;

import net.sf.samtools.*;


import java.io.*;

public class Main {

	public static void main(String[] args) {

		// ---------------------------------------------------------------------
		String fastaRef= "test_data/chr7.fa";
		String bam= "test_data/ds051.actb.bam";
		// ---------------------------------------------------------------------
		
		System.err.print("Reading fasta reference... ");
		FastaReference fastaFile= new FastaReference(new File(fastaRef)); 
		System.err.println("Done");
		
		SAMFileReader samfile= new SAMFileReader(new File(bam));
		
		int i= 0;
		ReadProfile readProfile= new ReadProfile();
		for (SAMRecord rec : samfile){
			AlignedRead alignedRead= new AlignedRead(rec, fastaFile);
			readProfile.addReadMethyl(alignedRead);
			i++;
			if (i % 1000000 == 0){
				System.err.println(i);
			}
		}
		System.err.println("N. reads: " + i);
		System.out.println(readProfile.toString());
		
		samfile.close();

	}
}