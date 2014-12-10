package sequenceMatcher;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;


public class ArgParse {
	
	public static String VERSION= "0.1.0";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser("SequenceMatcher")
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
						+ "This program answers the question: Which sequences in a file are "
						+ "similar to the sequences in another file. Where 'similarity' is the edit distance.\n"
						+ "\n"
						+ "Output format:\n"
						+ "1. Sequence name from the A file\n"
						+ "2. Sequence name from the B file\n"
						+ "3. Strand (+: Sequence matched as they are; -: Sequence A reverse complemented)\n"
						+ "4. Edit distance\n"
						+ "5. Length sequence A\n"
						+ "6. Length sequence B\n"
						+ "7. Sequence A, if matched as rev. comp. it is rev comp'd here.\n"
						+ "8. Sequence B\n");
		parser.addArgument("-a", "--a")
			.type(String.class)
			.required(true)
			.help("Fasta file 'A'; can be gzip'd. Use '-' to read from stdin");
		parser.addArgument("-b", "--b")
			.type(String.class)
			.required(true)
			.help("Fasta file 'B'; can be gzip'd. Use '-' to read from stdin");
		parser.addArgument("-m", "--method")
			.type(String.class)
			.required(false)
			.setDefault("Leven")
			.help("Method to determine the edit distance: Leven (Levenshtein, the default) or Hamming"
					+ "Humming is much faster.");
		parser.addArgument("-nm", "--nm")
			.type(Integer.class)
			.setDefault(-1)
			.help("Maximum edit distance to output a match. If -1 (default) all sequence pairs are returned.");
		parser.addArgument("-norc", "--norc")
			.action(Arguments.storeTrue())
			.help("Do not reverse complement the sequences in file A. I.e. only match sequences as they are.");	
		
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
		//		
	}
}
