/***************************************************************** 
*COMP 550 SI6 - Computing above/below for each flag 

Input : Red and blue segments with no red/red or blue/blue crossing 
	and no endpoints inside segments.
Output: For each flag, the segment above and below of the same and 
	the opposite color

Date: Oct 15th, 2014 by  Wenhua Guan

******************************************************************/
#include <stdio.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <limits>
#include <string>
#include <set>

using namespace std;

//Define 1st class: Point
class Point {

	public:	int x, y;
	
	Point(int x0, int y0){
		x = x0;
		y = y0;
	}
	Point(Point *P0){
		x = P0->x;
		y = P0->y;
	}
	
	~Point(){}

	//Overload operator =
	Point& operator=(const Point& P0){
		x = P0.x;
		y = P0.y;

		return *this;
	
	}
	
	//Comparison by x, breaking ties by y
	bool LessThan(Point *q){
		return (x < q->x)||((x == q->x) && (y < q->y));
	}
	//Define two points equal
	bool Equal(Point *q){
		return((x == q->x) && (y == q->y));
	}
	// To check whether P0P2 is on the left of P0P1
	static int Orientation(Point *P0, Point *P1, Point *P2){

	/*	long area = (P2->x - P0->x)*(P1->y - P0->y) - 
			    (P1->x - P0->x)*(P2->y - P0->y); */
		//Nov 30 change
                long area = (P1->x - P0->x)*(P2->y - P0->y) -
			    (P2->x - P0->x)*(P1->y - P0->y);

		if (area > 0){
			return 1; 	// left orientation
		}else if (area < 0){
			return -1; 	// right orientation
		}else{
			return 0; 	//colinear
		}
	}
};

class flag;
class Segment{
 
	//Members
	public:	flag  *start, *end;
	
	//constructor
	Segment(flag *start0, flag *end0)
	{
		start = start0;
		end = end0;  
	} 
	~Segment(){}

	//Comparision :P is on which Side of the seg
	static int Side(Segment *, Point *);

	//Check intersection
	static bool SEGEMENT_INTERSECT(Segment *, Segment *);

	static double Slope(flag *, flag *);

};

class flag {
	public: Point *P;
		Segment *Seg;
		int NoSeg;
		int typetag;
		string color;
			
	flag(Point *P0,int NoSeg0, int typetag0, string color0){
		P = P0;
		NoSeg = NoSeg0;
		typetag = typetag0;
		color = color0;
	}
	// set value to segment
	void SetSegment(Segment *Seg0){
		Seg = Seg0;	
	}
	~flag(){}
};

int Segment::Side(Segment *Seg, Point *P){
//	if(!P->Equal(Seg->start->P)&&(!P->Equal(Seg->end->P))){
	//Check P is above or belwo seg
		int side = Point :: Orientation(Seg->start->P, Seg->end->P, P);

		if(side > 0){
			return 1;  //Above Nov 30	// below 
		}else if( side < 0){
			return -1;  //Below Nov 30	// above
		}else{
			return 0;  	// colinear
		}
//	}

}
bool Segment::SEGEMENT_INTERSECT(Segment *Seg1, Segment *Seg2){
		int d1 = Side(Seg1, Seg2->start->P);
		int d2 = Side(Seg1, Seg2->end->P);
		int d3 = Side(Seg2, Seg1->start->P);
		int d4 = Side(Seg2, Seg1->end->P);

		if(d1*d2<=0 && d3*d4<=0){
			return true;
		}else{
			return false;
		}
	}	

double  Segment:: Slope(flag *start, flag *end){
			 if(start->P->x != end->P->x){
				double slope = (double)(end->P->y - start->P->y)/
					(double)(end->P->x - start->P->x);
				return slope;
			 }else{
				return numeric_limits<double>::max();	 
		 	 }
}

// Compare flags, return true if flag f0 < f1, vise versa
 bool compareflags(const flag *f0, const flag *f1){
	if(f0->P->Equal(f1->P)){
		cout<<f0->NoSeg<<f0->color<<f0->typetag<<":"<<f1->NoSeg<<f1->color<<f1->typetag<<
		"("<<Segment :: Slope(f0->Seg->start, f0->Seg->end)<<","<<
		Segment::Slope(f1->Seg->start, f1->Seg->end)<<")"<<endl;
	
	}

 	if (!(f0->P->Equal(f1->P))){
		return f0->P->LessThan(f1->P);// Point: Compare by x, breaking ties by y
	}else if (f0->typetag != f1->typetag){
		return f0->typetag < f1->typetag;            // Tyep : terminal < start
	}else if(Segment::Slope(f0->Seg->start, f0->Seg->end) !=
		 Segment::Slope(f1->Seg->start, f1->Seg->end) ){
		return (Segment :: Slope(f0->Seg->start, f0->Seg->end)<
		Segment::Slope(f1->Seg->start, f1->Seg->end)); // Slope
	}else if(f0->typetag == 1){ //start blue<  start red  
		return f0->color <  f1->color; 
	}else if(f0->typetag == 0){ // red term < blue term
		return f0->color > f1->color;
	}else{
		 cout <<"Overlapping lines of same color error!" <<endl;
		 return false;
	}
}

//EventLess is the operator of the tree structure Set
/*==============================================================================
   seg1 < seg2
   1.if seg1->start->x < seg2->start->x, Side(Seg1,Seg2->start)<0 => seg1 < seg2
   2.if seg1->start->x > seg2->start->x, Side(Seg1,Seg2->start)<0 => seg1 > seg2
   3 same x for start, then compare y value
   4 same start:  slope1 < slope 2 => seg1 < seg2 
==============================================================================*/   
struct EventLess{
	bool operator()(Segment *Seg1, Segment *Seg2){
		//If Seg1 < Seg2, return true; Vise versa.
		cout<<"compare"<<Seg1->start->NoSeg<<Seg1->start->color<<","<<
			Seg2->start->NoSeg<<Seg2->start->color<<endl;
		if( Seg1->start->P->x < Seg2->start->P->x ){
			cout<<"Seg1x<seg2x:"<<Segment::Side(Seg1,Seg2->start->P)<<endl;
			return (Segment::Side(Seg1,Seg2->start->P) <  0);
		}else if(Seg1->start->P->x > Seg2->start->P->x ){
			cout<<"Seg2x>Seg1x:"<<Segment::Side(Seg1,Seg2->start->P)<<endl;
			return (Segment::Side(Seg2,Seg1->start->P) > 0);
		}else if(Seg1->start->P->y != Seg2->start->P->y){
			cout<<"y:"<<Seg1->start->P->y<<"? "<<Seg2->start->P->y<<endl;
			return Seg1->start->P->y < Seg2->start->P->y;
		}else{
			double slope1 = Segment:: Slope(Seg1->start,Seg1->end);
			double slope2 = Segment:: Slope(Seg2->start,Seg2->end);			
			cout<<"Slope:"<<slope1<<"?"<<slope2<<endl;
			return slope1 <  slope2; // change to slope instead of -slope Nov. 29
		}
	}
};

// Insert new segment if there is end point flagtag->P in the middle of segtag.
void InsertNewSeg(vector<Segment *>& Seg, Segment* segtag, flag* flagtag){
	Point *Pstart = new Point(flagtag->P);
	Point *Pend = new Point(segtag->end->P);
				
	flag *sflag = new flag(Pstart, Seg.size(),1,segtag->start->color);
	flag *eflag = new flag(Pend, Seg.size(),0,segtag->end->color);
				
	Segment *s = new Segment(sflag,eflag);
	sflag->SetSegment(s);
	eflag->SetSegment(s);
				
	Seg.push_back(s);
			
	*(segtag->end->P) = *(flagtag->P);

}

int main(int argc, char *argv[]){
	int m, n, k;
	
	//Read coordinates of segments into variables
	if(argc !=3){
		cout<<"ERROR for running command!"<<endl;
		exit(EXIT_FAILURE);
	}
	
	FILE *fp;
	fp = fopen(argv[1],"r");
	if(fp ==NULL){
		cout<<"ERROR: Can't open input file!"<<endl;
		exit(EXIT_FAILURE);
	}
    //Input Iteration number
	int runtimes = atoi(argv[2]);
	int jointflag = 0;    // jointflag = 1 if there is joint segments 
	
	//Begin to read data
	fscanf(fp,"%d %d %d\n",&m, &n, &k );
      	
	vector<Segment *> RedSeg;
	vector<Segment *> BlueSeg;
	vector<flag *> FlagList;
	
	//For tight Sentinel line
	int xmin, ymin, xmax, ymax;
 	xmin = numeric_limits<int>::max();
	ymin = numeric_limits<int>::max();
	xmax = numeric_limits<int>::min();
	ymax = numeric_limits<int>::min(); // */
	
	//1. Read red segment
	for (int i = 0; i < m; i++){
		int x1,y1,x2,y2;
		fscanf(fp,"%d %d %d %d\n", &x1,&y1, &x2, &y2);
		
 		xmin = min(xmin, x1); xmin = min(xmin, x2);
		xmax = max(xmax, x1); xmax = max(xmax, x2);
		ymin = min(ymin, y1); ymin = min(ymin, y2);
		ymax = max(ymax, y1); ymax = max(ymax, y2);
		Point *P = new Point(x1, y1), *Q = new Point(x2, y2);
		//Nov30 start is always on the left of end
		if(x1>x2){
			P->x = x2;
			P->y = y2;
	                Q->x = x1;
	                Q->y = y1;
		}

		//start: typetag =1; end: typetag=0
		flag *start = new flag(P,i+1,1,"R"), *end = new flag(Q,i+1,0,"R");
		Segment *s = new Segment(start, end);
		start->SetSegment(s);		end->SetSegment(s);
		
		RedSeg.push_back(s); 		
		FlagList.push_back(start);   
		FlagList.push_back(end); 
	}
		
	// 2. Read blue segments
	for(int j = 0; j < n; j++){
		int x1,y1,x2,y2;
		fscanf(fp,"%d %d %d %d\n", &x1,&y1, &x2, &y2);	
		xmin = min(xmin, x1); xmin = min(xmin, x2);
		xmax = max(xmax, x1); xmax = max(xmax, x2);
		ymin = min(ymin, y1); ymin = min(ymin, y2);
		ymax = max(ymax, y1); ymax = max(ymax, y2);
		Point *P = new Point(x1,y1);
		Point *Q = new Point(x2,y2);
		// Nov 30 added exchange. start is always on the left of end
		if(x1>x2){
			P->x = x2;
			P->y = y2;
			Q->x = x1;
			Q->y = y1;
		}

		flag *start = new flag(P,j+1,1,"B"), *end = new flag(Q, j+1, 0,"B");  
	        Segment *s = new Segment(start, end);
		start->SetSegment(s); 	end->SetSegment(s);

		BlueSeg.push_back(s);
		FlagList.push_back(start);
		FlagList.push_back(end);

	}

	fclose(fp);
	// Tight Sentinel line for Red Segments: Red below and Red above 
	
 	Point RBS(xmin- 5, ymin-5), RBT(xmax+5, ymin-5);  // Red below start/terminal 
	Point RAS(xmin- 5, ymax+5), RAT(xmax+5, ymax+5);  // Red above start/terminal
	flag RBSF(&RBS, 0, 1, "R"), RBTF(&RBT, 0, 0, "R");  //Red Below start/terminal flag
	flag RASF(&RAS, m+1, 1, "R"), RATF(&RAT, m+1, 0, "R");   // Red above start/terminal flag
	Segment Rsegbotom(&RBSF,&RBTF), Rsegtop(&RASF,&RATF);   //Red bottom/top sentinel line */
	
	Point BBS(xmin-10, ymin-10), BBT(xmax+10, ymin-10);  // Blue below start/terminal 
	Point BAS(xmin-10, ymax+10), BAT(xmax+10, ymax+10);  // Blue above start/terminal
	flag BBSF(&BBS, 0, 1, "B"), BBTF(&BBT, 0, 0, "B");  //Blue Below start/terminal flag
	flag BASF(&BAS, n+1, 1, "B"), BATF(&BAT, n+1, 0, "B");   //Blue above start/terminal flag
	Segment Bsegbotom(&BBSF,&BBTF), Bsegtop(&BASF,&BATF);   //Blue bottom/top sentinel line
   
	cout<<"Rsegbotom: ("<<RBS.x<<","<<RBS.y<<")("<<RBT.x<<","<<RBT.y<<")"<<endl;
	cout<<"Rsegtop:    ("<<RAS.x<<","<<RAS.y<<")("<<RAT.x<<","<<RAT.y<<")"<<endl;
	cout<<"Bsegbotom: ("<<BBS.x<<","<<BBS.y<<")("<<BBT.x<<","<<BBT.y<<")"<<endl;
	cout<<"Bsegtop:    ("<<BAS.x<<","<<BAS.y<<")("<<BAT.x<<","<<BAT.y<<")"<<endl;

	// Sort FlagList
	sort(FlagList.begin(), FlagList.end(),compareflags);
	for (vector<int>:: size_type i =0; i < FlagList.size(); i++){
		cout << FlagList[i]->NoSeg<<FlagList[i]->color
		<<(FlagList[i]->typetag ==0 ? "T":"S")<<
		"("<<FlagList[i]->P->x<<","<<FlagList[i]->P->y<<")"<<endl;
  	}
	
	// Sweep algorithm
	set<Segment*, EventLess> Event;
	set<Segment*, EventLess> :: iterator Above, Below;
	
	//Insert the sentinel line
	Event.insert(&Rsegbotom);
	Event.insert(&Rsegtop);

	Below = Event.lower_bound(&Rsegtop);
	Below--;
	cout<<"Below of Rsegtop "<<(*Below)->start->NoSeg<<(*Below)->start->color<<"("<<
		(*Below)->start->P->x <<","<<(*Below)->start->P->y<<")"<<"("<<
		(*Below)->end->P->x <<","<<(*Below)->end->P->y<<")"<<endl;
	
	Event.insert(&Bsegbotom);
	
	Above = Event.upper_bound(&Bsegbotom);
	 cout<<"Above of Bsegbotom "<<(*Above)->start->NoSeg<<(*Above)->start->color<<"("<<
	      (*Above)->start->P->x <<","<<(*Above)->start->P->y<<")"<<"("<<
	      (*Above)->end->P->x <<","<<(*Above)->end->P->y<<")"<<endl;

	Event.insert(&Bsegtop);
	Below = Event.lower_bound(&Bsegtop);
	Below--;
	cout<<"Below of Bsegtop "<<(*Below)->start->NoSeg<<(*Below)->start->color<<"("<<
		(*Below)->start->P->x<<"," <<(*Below)->start->P->y<<")"<<"("<<
		(*Below)->end->P->x<<"," <<(*Below)->end->P->y<<")"<<endl;

	cout<<"flag"<<" "<<"sb"<<" "<<"sa"<<" "<<"ob"<<" "<<"oa"<<endl;
	
	//Run the code below for many times
	for(int k = 0; k < runtimes; k++){

        //Use sweep line algorithm to break segments 
		//Sweep the line through flags
		for(vector<int>::size_type i=0; i < FlagList.size(); ++i){
			int sb, sa, ob, oa;
			//start
			if(FlagList[i]->typetag == 1){
				Event.insert(FlagList[i]->Seg);
				Below = Event.lower_bound(FlagList[i]->Seg);
				Above = Event.upper_bound(FlagList[i]->Seg);
				Below--;

				// get sb and ob
				if((*Below)->start->color == FlagList[i]->color){
					sb = (*Below)->start->NoSeg; //same color
					while((*Below) && ((*Below)->start->color == FlagList[i]->color)){
						Below--;
					}
					ob = (*Below)->start->NoSeg;
				}else{
					ob = (*Below)->start->NoSeg; //opposite color
					while((*Below) && ((*Below)->start->color != FlagList[i]->color)){
						Below--;
					}
					sb = (*Below)->start->NoSeg;
				}				
				// get sa and oa
				if((*Above)->start->color == FlagList[i]->color){
					sa = (*Above)->start->NoSeg; //same color
					while( (*Above) && ((*Above)->start->color == FlagList[i]->color)){
						Above++;
					}
					oa = (*Above)->start->NoSeg;
				}else{
					oa = (*Above)->start->NoSeg; //opposite color
					while((*Above) && ((*Above)->start->color != FlagList[i]->color)){
						Above++;
					}
					sa = (*Above)->start->NoSeg;
				}
				
				cout<<FlagList[i]->NoSeg<<FlagList[i]->color<<
				(FlagList[i]->typetag ==0 ? "T":"S")<<"   "<<sb<<"  "<<sa<<"  "<<ob<<"  "
				<<oa<<endl;

			}

			//ternimal
			if(FlagList[i]->typetag==0){
				Below = Event.lower_bound(FlagList[i]->Seg);
				Above = Event.upper_bound(FlagList[i]->Seg);
				Below--;
			
				sb = FlagList[i]->NoSeg;
				sa = FlagList[i]->NoSeg;
				//Get ob and oa for terminal
				while( (*Below) && ((*Below)->start->color == FlagList[i]->color)){
						Below--;
				}
				ob = (*Below)->start->NoSeg;

				while((*Above) && ((*Above)->start->color == FlagList[i]->color)){
						Above++;
				}
				oa = (*Above)->start->NoSeg;
				
				cout<<FlagList[i]->NoSeg<<FlagList[i]->color<<
				(FlagList[i]->typetag ==0 ? "T":"S")<<"   "<<sb<<"  "<<sa<<"  "<<ob<<"  "
				<<oa<<endl;
				Event.erase(FlagList[i]->Seg);
			}
	
		}
		Event.clear();

	}

	
	//Delete elements in vector
	for(vector<int>::size_type i=0; i < RedSeg.size(); i++){ 
		if(RedSeg[i]){	
			delete RedSeg[i]->start->P;
			delete RedSeg[i]->end->P;
			delete RedSeg[i]->start;
			delete RedSeg[i]->end;
			delete RedSeg[i];
		}

	}

	for(vector<int>::size_type j=0; j < BlueSeg.size(); j++){
		if(BlueSeg[j]){
			delete BlueSeg[j]->start->P;
			delete BlueSeg[j]->end->P;
			delete BlueSeg[j]->start;
			delete BlueSeg[j]->end;
			delete BlueSeg[j];
		}
	}
	
    // Free memory
	vector<Segment *>().swap(RedSeg);
	vector<Segment *>().swap(BlueSeg);
	vector<flag *>().swap(FlagList);

	return 0;

}

	  
