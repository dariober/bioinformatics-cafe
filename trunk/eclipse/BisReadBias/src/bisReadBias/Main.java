package bisReadBias;

import net.sf.samtools.*;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.*;

import bisReadBias.ArgParse;

public class Main {

	public static void main(String[] args) {

		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		
		String fasta= opts.getString("fasta");
		String input= opts.getString("input");
		int step= opts.getInt("step");
		int stopAfter= opts.getInt("stopAfter");
		
		/* Adjust arguments */
		if (step <= 0){
			step= 1;
		} 
		
		/* ------------------------------------------------------------------- */
		System.err.print("Reading fasta reference... ");
		FastaReference fastaFile= new FastaReference(new File(fasta)); 
		System.err.println("Done");
		
		/* ------------------------------------------------------------------- */
		SAMFileReader samfile;
		if (input == "-"){
			samfile= new SAMFileReader(System.in);
		} else {
			samfile= new SAMFileReader(new File(input));
		}
		
		/* ------------------------------------------------------------------- */
		int nrec= 0;
		int ncollected= 0;
		ReadProfile readProfile= new ReadProfile();
		for (SAMRecord rec : samfile){
			if(nrec % step == 0){
				AlignedRead alignedRead= new AlignedRead(rec, fastaFile);
				readProfile.addReadMethyl(alignedRead);
				ncollected++;
			}
			nrec++;
			if(stopAfter > 0 && ncollected >= stopAfter){
				break;
			}
			if (nrec % 1000000 == 0){
				System.err.println(nrec);
			}
		}
		System.err.println("N. reads: " + nrec + "; N. collected: " + ncollected);
		System.out.println(readProfile.toString());
		
		samfile.close();

	}
}