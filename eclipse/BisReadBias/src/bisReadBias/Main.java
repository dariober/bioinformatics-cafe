package bisReadBias;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.*;

import bisReadBias.ArgParse;

public class Main {

	public static void main(String[] args) throws IOException, InterruptedException {

		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		
		String fasta= opts.getString("fasta");
		String input= opts.getString("input");
		String output= opts.getString("output");
		int step= opts.getInt("step");
		int stopAfter= opts.getInt("stopAfter");
		
		/* Adjust arguments */
		if (step <= 0){
			step= 1;
		} 
		if (output.equals("")){
			if (input.equals("-"))
				output= "bisReadBias";
			else {
				output= input.replaceAll("sam$|bam$", "") + "bisReadBias";
			}
		}
		
		/* ------------------------------------------------------------------- */
		System.err.print("Reading fasta reference... ");
		FastaReference fastaFile= new FastaReference(new File(fasta)); 
		System.err.println("Done");
		
		/* ------------------------------------------------------------------- */
		SAMFileReader samfile;
		if (input.equals("-")){
			samfile= new SAMFileReader(System.in);
		} else {
			samfile= new SAMFileReader(new File(input));
		}
		samfile.setValidationStringency(ValidationStringency.SILENT);
		
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
		
		FileWriter fw = new FileWriter(output + ".txt");
		fw.write(readProfile.toString());
		fw.close();
		
		ScriptRunner sr= new ScriptRunner(); 
		sr.runRscript(sr.getTmpFileFromResourceFile("Plotter.R"), output + ".txt");
		
		samfile.close();

	}
}