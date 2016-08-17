#BOOST=/nfs/wsi/bs/share/usr/feldhahn/Projects/immunopeptidomics/Trie/boost_1_39_0/
#BOOST=/share/opt/x86_64_sl6/boost_1_53_0/
#CFLAGS= -I/opt/local/include -O3 -I$(BOOST) -Wno-deprecated 
CFLAGS=-I/usr/include/boost/ -O3 -Wall
CXXFLAGS=$(CFLAGS)
#LDFLAGS=-lboost_serialization

LDFLAGS=-L/usr/lib/x86_64-linux-gnu/ -lboost_serialization
#LDFLAGS=-L/share/opt/x86_64_sl6/boost_1_53_0/lib/ -lboost_serialization
#LDFLAGS=-L/share/usr/feldhahn/Projects/immunopeptidomics/Trie/boost_1_39_0/lib/ -lboost_serialization-gcc43-mt
#LDFLAGS=-L/usr/lib64/ -lboost_serialization
#LDFLAGS=-L/share/usr/feldhahn/Projects/immunopeptidomics/Trie/boost_1_39_0/lib/ -lboost_serialization

all:	compute_distance get_TrieArray compute_distances_ivac debug generate_Peptides

compute_distance.o: compute_distance.C Matrix.h Sequences.h TrieArray.h Trie.h 

compute_distance:	compute_distance.o
	$(CXX) $(LDFLAGS) $< -o compute_distance

compute_distances_ivac.o: compute_distances_ivac.C Matrix.h Sequences.h TrieArray.h Trie.h 

compute_distances_ivac:	compute_distances_ivac.o
	$(CXX) $(LDFLAGS) -o compute_distances_ivac  $<


test_map.o: test_map.C Matrix.h Sequences.h TrieArray.h Trie.h 

test_map:	test_map.o
	$(CXX) $(LDFLAGS) -o test_map $<

get_TrieArray.o: get_TrieArray.C Matrix.h Sequences.h TrieArray.h Trie.h 	

get_TrieArray: get_TrieArray.o
	$(CXX)  $(LDFLAGS) -o get_TrieArray $<

debug.o: debug.C Matrix.h Sequences.h TrieArray.h Trie.h

debug: debug.o
	$(CXX)  $(LDFLAGS) -o debug $<
	
generate_Peptides.o: generate_Peptides.C Matrix.h Sequences.h TrieArray.h Trie.h

generate_Peptides: generate_Peptides.o
	$(CXX)  $(LDFLAGS) -o generate_Peptides $<
#----------------------------------------------------------------------------------


closest_peptide_pair.o:	closest_peptide_pair.C Matrix.h Sequences.h TrieArray.h Trie.h 	

closest_peptide_pair:	closest_peptide_pair.o
	$(CXX) $(LDFLAGS) -o closest_peptide_pair $<
	
closest_peptide_pair_proteome.o:	closest_peptide_pair_proteome.C Matrix.h Sequences.h TrieArray.h Trie.h 	

closest_peptide_pair_proteome:	closest_peptide_pair_proteome.o
	$(CXX) $(LDFLAGS) -o closest_peptide_pair_proteome $<
	
closest_peptide_complete.o:	closest_peptide_complete.C Matrix.h Sequences.h TrieArray.h Trie.h 	

closest_peptide_complete:	closest_peptide_complete.o
	$(CXX) $(LDFLAGS) -o closest_peptide_complete $<

closest_peptide_w1.o:	closest_peptide_w1.C Matrix.h Sequences.h TrieArray.h Trie.h 	

closest_peptide_w1:	closest_peptide_w1.o
	$(CXX) $(LDFLAGS) -o closest_peptide_w1 $<

closest_peptide_w2.o:	closest_peptide_w2.C Matrix.h Sequences.h TrieArray.h Trie.h 	

closest_peptide_w2:	closest_peptide_w2.o
	$(CXX) $(LDFLAGS) -o closest_peptide_w2 $<


#-------------------------------------------------------
closest_peptide_A2402.o:	closest_peptide_A2402.C Matrix.h Sequences.h TrieArray.h Trie.h 
closest_peptide_A2402:	closest_peptide_A2402.o
	$(CXX) $(LDFLAGS) -o closest_peptide_A2402 $<
	
closest_peptide_A2402_w1.o:	closest_peptide_A2402_w1.C Matrix.h Sequences.h TrieArray.h Trie.h 
closest_peptide_A2402_w1:	closest_peptide_A2402_w1.o
	$(CXX) $(LDFLAGS) -o closest_peptide_A2402_w1 $<
	
closest_peptide_A2402_w2.o:	closest_peptide_A2402_w2.C Matrix.h Sequences.h TrieArray.h Trie.h 
closest_peptide_A2402_w2:	closest_peptide_A2402_w2.o
	$(CXX) $(LDFLAGS) -o closest_peptide_A2402_w2 $<
	

	
#-------------------------------------------------------






sanity_check.o:	sanity_check.C Matrix.h Sequences.h TrieArray.h Trie.h TrieFunctions.h

sanity_check:	sanity_check.o
	$(CXX) $(LDFLAGS) -o sanity_check $<


check_pair.o: check_pair.C Matrix.h Sequences.h TrieArray.h Trie.h TrieFunctions.h

check_pair: check_pair.o
	$(CXX) $(LDFLAGS) -o check_pair $<


#demo: 
	#$(CXX) $(CXXFLAGS) -c demo.C -o demo.o
	#$(CXX) $(LDFLAGS) -o demo demo.o

