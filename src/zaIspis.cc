//za ispis
void opcenitInfo(std::vector<minimizer>& start_minimizers, std::vector<minimizer>& end_minimizers,
  std::vector<minimizer_hit_t>& start_hits_same, std::vector<minimizer_hit_t>& start_hits_rev,
  std::vector<minimizer_hit_t>& end_hits_same, std::vector<minimizer_hit_t>& end_hits_rev, unsigned int fobj, unsigned int sequence_length){

    printf("Br. %d\t", fobj);
    printf("velicina sekvence: %5u\t", sequence_length);
    printf("Velicina: start_minimizers %lu\tend_minimizers %lu\t", start_minimizers.size(), end_minimizers.size());
    printf("Hitova: start_same %lu\tstart_rev %lu\tend_same %lu\tend_rev %lu\n", start_hits_same.size(), start_hits_rev.size(),
          end_hits_same.size(), end_hits_rev.size());
}

void regijeInfo(std::unordered_map<unsigned int, int>& start_hits_region, std::unordered_map<unsigned int, int>& start_hits_region_rev, 
  std::unordered_map<unsigned int, int>& end_hits_region, std::unordered_map<unsigned int, int>& end_hits_region_rev){

    printf("Pocetak:\n");
    for (auto& it: start_hits_region) {
      if(it.second < 5) continue;
      printf("Regija %d ima %d hitova\n", it.first, it.second);
    }
    printf("Reverse------------------\n");
    for (auto& it: start_hits_region_rev) {
      if(it.second < 5) continue;
      printf("Regija %d ima %d hitova\n", it.first, it.second);
    }
    printf("Kraj: \n");
    for (auto& it: end_hits_region) {
      if(it.second < 5) continue;
      printf("Regija %d ima %d hitova\n", it.first, it.second);
    }
    printf("Reverse------------------\n");
    for (auto& it: end_hits_region_rev) {
      if(it.second < 5) continue;
      printf("Regija %d ima %d hitova\n", it.first, it.second);
    }
}

void top3Info(std::vector<region_hits>& start_hits_top, std::vector<region_hits>& start_hits_top_rev,
  std::vector<region_hits>& end_hits_top, std::vector<region_hits>& end_hits_top_rev){
    printf("Pocetak:\n");
    for(auto it : start_hits_top){
      printf("Regija %d ima %d hitova\n", it.first, it.second);
    }
    for(auto it : start_hits_top_rev){
      printf("Regija %d ima %d hitova, reverse\n", it.first, it.second);
    }
    printf("Kraj:\n");
    for(auto it : end_hits_top){
      printf("Regija %d ima %d hitova\n", it.first, it.second);
    }
    for(auto it : end_hits_top_rev){
      printf("Regija %d ima %d hitova, reverse\n", it.first, it.second);
    }
}

    // printf("%s\t%u\t%u\t%u\t%c\t%u\t%u\n",fastaq_objects1[fobj]->name.c_str(), sequence_length, start, end, rev ? '-' : '+', ref_start, ref_end);

    // ispis ako se zeli vidjet podatci za regije i kak se odreduje
    // printf("%s %lu \n", fastaq_objects1[fobj]->name.c_str(), sequence_length);
    // regijeInfo(start_hits_region, start_hits_region_rev, end_hits_region, end_hits_region_rev);
    // printf("Top 3 uz region_size:%d i nade starting_region %lu, ending_region %lu\n", region_size, starting_region, ending_region);
    // top3Info(start_hits_top, start_hits_top_rev, end_hits_top, end_hits_top_rev);
    // printf("%u\t%u\t%u\t%c\t%u\t%u\n", sequence_length, start, end, rev ? '-' : '+', ref_start, ref_end);
    // printf("\n");