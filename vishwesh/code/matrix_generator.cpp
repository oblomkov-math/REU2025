#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>        // for std::strcmp
#include <cmath>
#include <regex>
#include <unordered_set>
#include <algorithm>
#include <cstdlib>
#include <limits>

#include "prime_power.hpp"

static const char* MATRIX_DIR = "../matrices/";

// Display a simple progress bar to stdout
void showProgress(std::uint64_t current, std::uint64_t total) {
    const int width = 50;
    float fraction = total ? static_cast<float>(current) / total : 1.0f;
    int filled = static_cast<int>(fraction * width);
    std::cout << "[";
    for (int i = 0; i < width; ++i) {
        if (i < filled) std::cout << '=';
        else if (i == filled) std::cout << '>';
        else std::cout << ' ';
    }
    std::cout << "] " << int(fraction * 100) << "%\r";
    std::cout.flush();
    if (current >= total) std::cout << std::endl;
}

// Convert a numeric ID into an N×N matrix over F_p in row-major order
// Most-significant digit first: for ID=1 and p=2,N=2 -> <0,0,0,1>
void idToMatrix(std::uint64_t id, std::uint64_t p, int N, std::vector<std::uint64_t>& out) {
    out.assign(N * N, 0);
    for (int i = 0; i < N * N; ++i) {
        out[N*N - 1 - i] = id % p;
        id /= p;
    }
}

// Multiply C = A * B over F_p, size N×N
void multiply(const std::vector<std::uint64_t>& A,
              const std::vector<std::uint64_t>& B,
              std::vector<std::uint64_t>& C,
              std::uint64_t p, int N) {
    std::fill(C.begin(), C.end(), 0);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < N; ++k) {
            auto aik = A[i*N + k];
            if (!aik) continue;
            for (int j = 0; j < N; ++j) {
                C[i*N + j] = (C[i*N + j] + aik * B[k*N + j]) % p;
            }
        }
    }
}

// Check nilpotency: M^N == 0
bool isNilpotent(const std::vector<std::uint64_t>& M, std::uint64_t p, int N) {
    std::vector<std::uint64_t> cur = M, tmp(N*N);
    for (int exp = 1; exp < N; ++exp) {
        multiply(cur, M, tmp, p, N);
        cur.swap(tmp);
        if (std::all_of(cur.begin(), cur.end(), [](auto v){ return v == 0; })) return true;
    }
    return false;
}

// Open or resume a matrix list (general/nilpotent)
void openMatrixList(const std::string& path,
                    std::fstream& fs,
                    bool& completed,
                    std::uint64_t& written,
                    std::uint64_t& nextId) {
    std::regex entryRe(R"(^\s*(\d+):)");
    completed = false;
    written = 0;
    nextId = 1;
    if (!std::ifstream(path)) {
        fs.open(path, std::ios::out);
        fs << "COMPLETE: 0\nCOUNT: -1\n";
        fs.flush();
    } else {
        fs.open(path, std::ios::in | std::ios::out);
        std::string line;
        std::getline(fs, line);
        completed = (line.find('1') != std::string::npos);
        std::getline(fs, line); // skip COUNT
        auto dataPos = fs.tellg();
        while (std::getline(fs, line)) {
            std::smatch m;
            if (std::regex_search(line, m, entryRe)) {
                nextId = std::stoull(m[1]) + 1;
                ++written;
            }
        }
        fs.clear();
        fs.seekp(dataPos);
        if (!completed && written > 0) std::cout << "Continuing from ID=" << nextId << " in " << path << std::endl;
    }
}

int main(int argc, char* argv[]) {
    bool gen=false, nil=false, commute=false;
    int N=0; std::uint64_t P=0;
    for (int i=1;i<argc;++i) {
        if (!std::strcmp(argv[i],"--general")) gen=true;
        else if (!std::strcmp(argv[i],"--nilpotent")) nil=true;
        else if (!std::strcmp(argv[i],"-c")) commute=true;
        else if (!N) N=std::atoi(argv[i]);
        else if (!P) P=std::strtoull(argv[i],nullptr,10);
        else {std::cerr<<"Unknown arg: "<<argv[i]<<std::endl;return 1;}
    }
    if ((gen&&nil)||(!gen&&!nil)) {std::cerr<<"Specify --general or --nilpotent."<<std::endl;return 1;}
    if (N<=0||P<2||!isPrimePower(P)) {std::cerr<<"Invalid N or P."<<std::endl;return 1;}

    std::string base=MATRIX_DIR;
    std::string fGen=base+"GENERAL_n="+std::to_string(N)+",p="+std::to_string(P)+".txt";
    std::string fNil=base+"NILPOTENT_n="+std::to_string(N)+",p="+std::to_string(P)+".txt";
    std::string listF=gen?fGen:fNil;
    std::string fCom=base+"COMMUTING_"+(gen?"GENERAL":"NILPOTENT")+"_n="+std::to_string(N)+",p="+std::to_string(P)+".txt";

    // auto-gen general for nilpotent or commute
    if ((nil||commute) && !gen) {
        std::ifstream chk(fGen);
        bool need=false;
        if (!chk) {std::cout<<"GENERAL missing; generating..."<<std::endl;need=true;}
        else {std::string l;std::getline(chk,l);if(l.find('1')==std::string::npos){std::cout<<"GENERAL incomplete; regenerating..."<<std::endl;need=true;}}
        if (need) {std::string cmd=std::string(argv[0])+" --general "+std::to_string(N)+" "+std::to_string(P);
            if (std::system(cmd.c_str())){std::cerr<<"Failed GENERAL"<<std::endl;return 1;}}
    }

    // generate/resume list
    std::fstream fsL; bool done=false; std::uint64_t cnt=0,nxt=1;
    openMatrixList(listF,fsL,done,cnt,nxt);
    if (!done) {
        __int128 t128=1;for(int i=0;i<N*N;++i) t128*=P;
        if (t128>std::numeric_limits<std::uint64_t>::max()){std::cerr<<"Too many matrices."<<std::endl;return 1;}
        auto total=(std::uint64_t)t128;
        std::vector<std::uint64_t> mat(N*N),tmp(N*N);
        for(std::uint64_t id=nxt;id<=total;++id){showProgress(id,total);idToMatrix(id,P,N,mat);
            if(gen||isNilpotent(mat,P,N)){
                fsL<<id<<": <";
                for(size_t k=0;k<mat.size();++k)fsL<<mat[k]<<(k+1<mat.size()?",":"");
                fsL<<">\n";++cnt;
            }
        }
        fsL.close();
        // headers COMPLETE then COUNT
        fsL.open(listF,std::ios::in|std::ios::out);
        fsL.seekp(0);
        fsL<<"COMPLETE: 1\n";
        fsL<<"COUNT: "<<cnt<<"\n";
        fsL.close();
        std::cout<<(gen?"GENERAL":"NILPOTENT")<<" matrix generation complete n="<<N<<", p="<<P<<"!"<<std::endl;
        done=true;
    } else {
        std::cout<<(gen?"GENERAL":"NILPOTENT")<<" matrix generation already complete n="<<N<<", p="<<P<<"!"<<std::endl;
    }

    if (commute) {
        if(!done){std::cerr<<"List incomplete; cannot commute."<<std::endl;return 1;}
        std::ifstream ifs(listF);
        std::string l;getline(ifs,l);getline(ifs,l);
        std::vector<std::uint64_t> ids;
        std::regex re(R"(^\s*(\d+):)");
        while(getline(ifs,l)){std::smatch m;if(std::regex_search(l,m,re))ids.push_back(std::stoull(m[1]));}
        ifs.close(); size_t M=ids.size();
        std::fstream fsC; bool cpl=false; std::uint64_t cp=0;
        std::unordered_set<std::uint64_t> seen;
        if(!std::ifstream(fCom)){fsC.open(fCom,std::ios::out);fsC<<"COMPLETE: 0\nCOUNT: -1\n";}
        else{fsC.open(fCom,std::ios::in|std::ios::out);getline(fsC,l);cpl=(l.find('1')!=std::string::npos);getline(fsC,l);
            auto pos=fsC.tellg(); std::regex pr(R"(^\s*\((\d+),(\d+)\))");
            while(getline(fsC,l)){std::smatch m;if(std::regex_search(l,m,pr)){auto a=std::stoull(m[1])-1,b=std::stoull(m[2])-1;seen.insert(a*M+b);++cp;}}
            fsC.clear();fsC.seekp(pos); if(cpl){std::cout<<(gen?"GENERAL":"NILPOTENT")<<" COMMUTING PAIR generation already complete n="<<N<<", p="<<P<<"!"<<std::endl;return 0;} 
            std::cout<<"Continuing commuting from idx="<<cp<<std::endl;
        }
        std::vector<std::uint64_t> A(N*N),B(N*N),AB(N*N),BA(N*N);
        auto totalP=M*M;
        for(size_t i=0;i<M;++i){idToMatrix(ids[i],P,N,A);
            for(size_t j=0;j<M;++j){auto idx=i*M+j; if(seen.count(idx)) continue; showProgress(idx,totalP);
                idToMatrix(ids[j],P,N,B);multiply(A,B,AB,P,N);multiply(B,A,BA,P,N);
                if(AB==BA){fsC<<"("<<ids[i]<<","<<ids[j]<<")\n";++cp;}seen.insert(idx);
            }}
        fsC.seekp(0);fsC<<"COMPLETE: 1\n"<<"COUNT: "<<cp<<"\n";fsC.close();
        std::cout<<(gen?"GENERAL":"NILPOTENT")<<" COMMUTING PAIR generation complete n="<<N<<", p="<<P<<"!"<<std::endl;
    }
    return 0;
}

