#include "visibility.h"

visibility::visibility(const vector<point>& pts)
{
	poly.resize(pts.size());
	for(int i = 0; i < pts.size(); i++)
		poly[i].set_point(pts[i].get_point_cartesian().first, 
						pts[i].get_point_cartesian().second);
}

vector<point> visibility::find_visibility_polygon(point cur_point)
{
	ret.clear();
	pts.clear();
	alpha.clear();
	pt = cur_point.get_point_cartesian();
	int chk = 0;
	flag = 0;
	int minpt = -1;
	long double dist = MAXVAL;

	for(int i = 0; i < poly.size(); i++)
	{
		long double val = twice_triangle_signed_area(poly[(i-1)%poly.size()], poly[i], cur_point);
		if(fabs(val) < eps)
		{
			chk = i;
			break;
		}
	}

	if(cur_point == poly[chk]) chk = (chk+1)%poly.size();
	
	if(chk)
	{
		if(!(cur_point == poly[chk] or cur_point == poly[chk-1]))
		{
			pts.push_back(point().get_point_cartesian());
			alpha.push_back(0.0);
		}
		else flag++;
	}
	else
	{
		//Find order
		for(int i = 0; i < poly.size(); i++)
		{
			long double tmp = (pt.first - poly[i].get_point_cartesian().first) *
							  (pt.first - poly[i].get_point_cartesian().first) +
							  (pt.second - poly[i].get_point_cartesian().second) *
							  (pt.second - poly[i].get_point_cartesian().second);
			if(tmp < dist)
			{
				dist = tmp;
				minpt = i;
			}
		}
		chk = minpt;
	}

	int idx = chk, cnt = 0;
	long double rad = point(poly[chk].get_point_cartesian().first-pt.first, 
							poly[chk].get_point_cartesian().second-pt.second).
							get_point_polar().second.get_value(true), prv = 0.0;
	
	do
	{
		auto tmp = point(poly[idx].get_point_cartesian().first-pt.first, 
						 poly[idx].get_point_cartesian().second-pt.second).get_point_polar();
		double val = tmp.second.get_value(true)-rad - prv;
		val += 4.0 * pi();
		
		while(val > pi()) val -= 2.0*pi();
		
		if(idx == chk and !cnt) 
		{
			alpha.push_back(0.0);
		}
		else
		{
			double sgn = 1.0;
			alpha.push_back(alpha[alpha.size()-1] + sgn*val);
			prv = tmp.second.get_value(true)-rad; 
		}
		
		tmp.second.set_angle(tmp.second.get_value(true)-rad);
		pts.push_back(point(tmp.first, tmp.second).get_point_cartesian());
		if(fabs(pts[pts.size()-1].first) < eps) pts[pts.size()-1].first = 0.0;
		if(fabs(pts[pts.size()-1].second) < eps) pts[pts.size()-1].second = 0.0;
		idx = (idx+1) % poly.size();
		
		if(idx == chk and minpt != -1) {chk++; minpt = -1; cnt++;}

	}while(idx != chk);

	if(cnt == 1)
	{
		if(alpha[alpha.size()-1] < 1.0) return ret;
		alpha.pop_back();
		alpha.push_back(2*pi());
	}
	
	int cur = 0, ccw;
	string op;
	
	while(!s.empty()) s.pop();
	s.push(0);
	//Nodes are visible hence advance 
	if(alpha[1] - alpha[0] > -eps) op = "advance";
	else 
	{
		op = "scan";
		ccw = 1;
		inf.set_point(MAXVAL, angle(0.0));
	}

	while(true)
	{
		if(op == "finish") break;
		else if(op == "advance") op = advance(cur, ccw);
		else if(op == "scan") op = scan(cur, ccw);
		else if(op == "retard") op = retard(cur, ccw);
	}

	while(!s.empty())
	{
		point tmpt = point(pts[s.top()].first,
						   pts[s.top()].second);
		double r = tmpt.get_point_polar().first;
		double ang = tmpt.get_point_polar().second.get_value(true);
		ang += rad;
		tmpt.set_point(r, angle(ang));
		
		double x = (fabs(tmpt.get_point_cartesian().first + pt.first) > eps) ?
						(tmpt.get_point_cartesian().first + pt.first): 0.0;
		double y = (fabs(tmpt.get_point_cartesian().second + pt.second) > eps) ?
						(tmpt.get_point_cartesian().second + pt.second): 0.0;				
		
		if(ret.empty() or point(x,y) != ret[ret.size()-1])
			ret.push_back(point(x, y));
		s.pop();

	}
	
	vector<point> ret_vec;
	
	for(int i = ret.size()-1; i >= 0; i--)
		ret_vec.push_back(ret[i]); 
	
	return ret_vec;
}

string visibility::advance(int& cur, int& ccw)
{
	//Handling of nodes in the advance face.
	while(true)
	{
		if((alpha[cur+1] - 2.0*pi()) < eps)
		{
			cur++;
			s.push(cur);

			if(cur + flag == poly.size()) return "finish";

			else if(alpha[cur+1] - alpha[cur] < -eps)
			{

				point A = point(pts[cur-1].first, pts[cur-1].second);
				point B = point(pts[cur].first, pts[cur].second);
				point C = point(pts[cur+1].first, pts[cur+1].second);
				long double turn = twice_triangle_signed_area(A, B, C);

				if(turn < 0)
				{
					ccw = true;
					inf.set_point(MAXVAL, angle(alpha[cur]));
					return "scan";
				}
				else
				{
					return "retard";
				}
			}
		}
		else
		{
			if(alpha[s.top()] - 2.0*pi() < -eps)
			{
				if(fabs(pts[cur].first-pts[cur+1].first) < eps)
				{
					pts.push_back(make_pair(pts[cur].first, 0.0));
					double val = 2.0*pi();
					while(val > alpha[cur+1]) val -= 2.0*pi();
					while(val < alpha[cur]) val += 2.0*pi(); 
					alpha.push_back(val);
					s.push(alpha.size()-1);
				}
				else
				{
					long double slp = (pts[cur].second-pts[cur+1].second);
					slp /= (pts[cur].first-pts[cur+1].first);
					long double c = pts[cur].second - slp*pts[cur].first;
					pts.push_back(make_pair(((-1.0*c)/slp), 0.0));
					alpha.push_back(2*pi());
					s.push(alpha.size()-1);
				}

				ccw = false;
				inf = point(pts[0].first, pts[0].second);
				return "scan";
			}
		}
	}
}

string visibility::retard(int& cur, int& ccw)
{
	//Handling of nodes in the retard face.
	while(true)
	{
		int top;
		while(true)
		{
			top = s.top();
			s.pop();
			
			if(alpha[s.top()] - alpha[cur+1] < -eps and (alpha[cur+1] - alpha[top]) < eps)
			{
				break;
			}
			long double tmp1 = twice_triangle_signed_area(point(pts[top].first, pts[top].second), 
												  point(pts[s.top()].first, pts[s.top()].second),
												  point(pts[cur+1].first, pts[cur+1].second));
			tmp1 *= twice_triangle_signed_area(point(pts[top].first, pts[top].second), 
				 					   point(pts[s.top()].first, pts[s.top()].second),
									   point(pts[cur].first, pts[cur].second));
			
			if(tmp1 > eps) continue;
			
			tmp1 *= twice_triangle_signed_area(point(pts[top].first, pts[top].second), 
				 					   	   point(pts[cur+1].first, pts[cur+1].second),
									       point(pts[cur].first, pts[cur].second));
			tmp1 *= twice_triangle_signed_area(point(pts[s.top()].first, pts[s.top()].second), 
				 					   	   		   point(pts[cur+1].first, pts[cur+1].second),
									               point(pts[cur].first, pts[cur].second));
			
			if(tmp1 < -eps) continue;
			
			if(fabs(alpha[top]-alpha[s.top()]) < eps and (alpha[cur+1] - alpha[s.top()]) < eps)
			{
				break;
			}
		}

		if(alpha[s.top()] - alpha[cur+1] < -eps)
		{
			cur++;
			point intr = intersection(point(), point(pts[cur].first, pts[cur].second),
									  		   point(pts[top].first, pts[top].second),
									  	point(pts[s.top()].first, pts[s.top()].second));
			
			double ang = intr.get_point_polar().second.get_value(true);
			while(ang > alpha[top]) ang -= 2.0*pi();
			while(ang < alpha[s.top()]) ang += 2.0*pi();
			
			pts.push_back(intr.get_point_cartesian());
			alpha.push_back(ang);
			
			s.push(alpha.size()-1);
			s.push(cur);
			
			point A = point(pts[cur-1].first, pts[cur-1].second);
			point B = point(pts[cur].first, pts[cur].second);
			point C = point(pts[cur+1].first, pts[cur+1].second);
			long double turn = twice_triangle_signed_area(A, B, C);
			
			if(cur + flag == poly.size()) return "finish";
			else if(alpha[cur+1] - alpha[cur] > -eps and turn < 0.0) return "advance";
			else if(alpha[cur+1] - alpha[cur] > eps and turn > 0.0)
			{
				ccw = 0;
				inf.set_point(pts[cur].first, pts[cur].second);
				s.pop();
				return "scan";
			}
			else s.pop();
		}
		else
		{
			point A = point(pts[cur-1].first, pts[cur-1].second);
			point B = point(pts[cur].first, pts[cur].second);
			point C = point(pts[cur+1].first, pts[cur+1].second);
			long double turn = twice_triangle_signed_area(A, B, C);
			
			if(fabs(alpha[s.top()]-alpha[cur+1]) < eps and alpha[cur+2] - alpha[cur+1] > eps and turn < 0.0)
			{
				cur++;
				s.push(cur);
				return "advance";
			}
			else
			{
				ccw =true;
				auto intr = intersection(point(pts[cur].first, pts[cur].second),
										point(pts[cur+1].first, pts[cur+1].second),
									  		point(pts[top].first, pts[top].second),
									point(pts[s.top()].first, pts[s.top()].second)).get_point_cartesian();
				inf.set_point(intr.first, intr.second);
				return "scan";
			}
		}
	}
}

string visibility::scan(int& cur, int& ccw)
{
	//Handling of nodes in the scan face.
	while(true)
	{
		cur++;
		long double tmp1 = twice_triangle_signed_area(inf, 
				point(pts[s.top()].first, pts[s.top()].second),
				    point(pts[cur+1].first, pts[cur+1].second));

		tmp1 *= twice_triangle_signed_area(inf, 
				point(pts[s.top()].first, pts[s.top()].second),
				    point(pts[cur].first, pts[cur].second));
		
		if(tmp1 > eps) continue;

		tmp1 *= twice_triangle_signed_area(inf, 
			 	point(pts[cur+1].first, pts[cur+1].second),
				    point(pts[cur].first, pts[cur].second));

		tmp1 *= twice_triangle_signed_area(point(pts[s.top()].first, pts[s.top()].second), 
			 					   	   		   point(pts[cur+1].first, pts[cur+1].second),
								               point(pts[cur].first, pts[cur].second));

		if(tmp1 < -eps) continue;

		if(ccw and alpha[cur+1] > alpha[s.top()] and alpha[s.top()] - alpha[cur] > -eps)
		{
			point intr = intersection(point(pts[cur+1].first, pts[cur+1].second),
				    				  point(pts[cur].first, pts[cur].second),
				    				  point(pts[s.top()].first, pts[s.top()].second),
				    				  inf);
			pts.push_back(intr.get_point_cartesian());
			
			long double ang = intr.get_point_polar().second.get_value(true);
			while(ang > alpha[cur+1]) ang -= 2.0*pi();
			while(ang < alpha[cur]) ang += 2.0*pi();
			alpha.push_back(ang);
			
			s.push(alpha.size()-1);
			return "advance";
		}
		if((!ccw) and alpha[cur+1]-alpha[s.top()] < eps and alpha[s.top()] < alpha[cur])
		{
			return "retard";
		}
	}
}