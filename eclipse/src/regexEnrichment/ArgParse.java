package regexEnrichment;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import java.util.*;

public class ArgParse {
	
	public static String VERSION= "0.1.0";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser("RegexEnrichment")
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
						+ "Split fasta file in windows to detect which windows "
						+ "are enriched for a regex pattern. "
						+ "I.e. where the pattern occurs more often than expected by random chance given "
						+ "the sequence composition of the tested window."
						+ "\n\n"
						+ "EXAMPLE:\n"
						+ "Scan chr1 for the forkhead motif in 1kb windows, sliding by 500bp. "
						+ "Read randomization results from file and save to the same "
						+ "file at the end of the run.\n\n"
						+ "java -Xmx2G -jar RegexEnrichment.jar -f chr1.hg19.fa.gz -p [GA]TAAA[CT]AA -w 1000 -s 500 -wr fox.ser -rt fox.ser"
						+ "\n\n"
						+ "OUTPUT:\n"
						+ "Bed file with one record per match (cols 1-3); 4th columns: record name; "
						+ "5th column: pvalue for the null hypothesis that the number of observed matches "
						+ "in the tested window comes from a random distribution given "
						+ "the observed nucleotide frequencies; "
						+ "6th: Strand of the match; 7-8th: coordinates of the tested window; "
						+ "9th: matched pattern, if --seqcase U or L, converted to upper or lower case, "
						+ "respectively; 10th (possibly deprecable): Number of hits in the randomized "
						+ "sequences."
						+ "");	
		parser.addArgument("-f", "--fasta")
			.type(String.class)
			.required(true)
			.help("Input fasta file to scan. Can be gzipped.");
		parser.addArgument("-o", "--output")
			.type(String.class)
			.required(false)
			.setDefault("-")
			.help("Output file. Use - to stdout");
		parser.addArgument("-p", "--pattern")
			.type(String.class)
			.required(false)
			.setDefault("([gG]{3,}\\w{1,7}){3,}[gG]{3,}")
			.help("Regular expression pattern to scan. Default is G-quadruplex motif");
		parser.addArgument("--casesensitive", "-cs")
			.action(Arguments.storeTrue())
			.help("Enables case sensitive matching. By default matching is "
					+ "case insensitive. See also --seqcase");
		
		parser.addArgument("-omc", "--orderMarkovChain")
			.type(Integer.class)
			.required(false)
			.setDefault(0)
			.help("Order of the Markov chain. Used to generate random sequences"
					+ " given the input sequence. Set to zero to disable the use "
					+ "of the Markov model and thus generate random sequences by pure "
					+ "random sampling with replacement");
		
		parser.addArgument("-pp", "--pctPseudocounts")
		.type(Double.class)
		.required(false)
		.setDefault(0.1)
		.help("Ignored if --orderMarkovChain=0: Assign to unobserved kmers in "
				+ "the Markov model pseudocounts to make them up this percentage of "
				+ "the total counts.");
		
		parser.addArgument("-w", "--window_size")
			.type(Integer.class)
			.required(false)
			.setDefault(1000)
			.help("Size of the sliding window.");
		parser.addArgument("--step", "-s")
			.type(Integer.class)
			.required(false)
			.setDefault(-1)
			.help("Step size for sliding window. Default (-1) sets to window size, i.e. "
					+ "windows are next to each other without overlap");		
		parser.addArgument("-r", "--nrand")
			.type(Integer.class)
			.required(false)
			.setDefault(100000)
			.help("Number of randomizations applied to each tested sequence window. "
					+ "The smallest non-zero pvalue will be 1/nrand.");
		parser.addArgument("-P", "--precision")
			.type(Integer.class)
			.required(false)
			.setDefault(20)
			.help("Precision to separate nucleotide frequencies. "
					+ "The randomization results for each encountered combination of "
					+ "nucleotide frequencies "
					+ "are kept in memory and they not re-computed if following windows "
					+ "have the same composition. "
					+ "--precision determines the degree with which similar compositions "
					+ "are treated as equivalent, therefore allowing the use of the same "
					+ "randomization results. "
					+ "Default 20 (1/20=0.05) means that frequencies are rounded to "
					+ "the nearest 1/20th. E.g. G:0.21 and G:0.23 are rounded to 0.2");
		parser.addArgument("--norevcomp")
			.action(Arguments.storeTrue())
			.help("Do not scan the reverse complement of the input sequence");
		String[] s= {"L", "U", "k"};
		parser.addArgument("--seqcase", "-c")
			.choices(Arrays.asList(s))
			.setDefault("k")
			.help("Convert the reference sequence to upper case (U), lower case (L) "
					+ "or do not change case (k). Converting to upper or lower case "
					+ "reduces the combinations of nucleotide frequencies thus speeding "
					+ "up the computation."
					+ "Forced to U if --casesensitive is not enabled otherwise default to k.");		
		parser.addArgument("--rLookUpTable", "-rt")
			.help("Read in the look up table of nucleotide frequencies and pattern matches "
					+ "from this file. If you are planning to scan different genomes for the same pattern "
					+ "and the same parameter arguments, save the randomization "
					+ "results to speed up the scanning.");		
		parser.addArgument("--wLookUpTable", "-wt")
			.help("Write out the look up table of nucleotide frequencies pattern matches "
					+ "to this file.");
		parser.addArgument("--verbose", "-V")
			.action(Arguments.storeTrue())
			.help("Print some progress report, once per minute.");
		parser.addArgument("--version", "-v").action(Arguments.version());
		Namespace opts= null;
		try{
			opts= parser.parseArgs(args);
		}
		catch(ArgumentParserException e) {
			parser.handleError(e);
			System.exit(1);
		}		
		return(opts);
	}
	
	public static void validateArgs(Namespace opts){
		
		if (opts.getInt("window_size") < 0 || 
			opts.getInt("nrand") < 0){
			
			System.out.println("Invalid args");
			System.exit(1);
		}
		if (opts.getInt("orderMarkovChain") < 0){
			System.out.println("Invalid orderMarkovModel: must be >= 0");
			System.exit(1);			
		}
		if (opts.getDouble("pctPseudocounts") < 0 || opts.getDouble("pctPseudocounts") > 1){
			System.out.println("Invalid pctPsedocounts: must be 0 < x < 1");
			System.exit(1);			
		}
	}
}
