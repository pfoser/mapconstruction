/*
Path-based graph distance 1.0
Copyright 2014 Mahmuda Ahmed, K. S. Hickmann and Carola Wenk

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

------------------------------------------------------------------------

This software is based on the following article. Please cite this
article when using this code as part of a research publication:

M. Ahmed, K. S. Hickmann, and C. Wenk.
Path-based distance for street map comparison.
arXiv:1309.6131, 2013.
------------------------------------------------------------------------

Author: Mahmuda Ahmed
Filename: IntervalComparator.java
 */
package mapmatchingbasics;
import java.util.Comparator;


public class IntervalComparator implements Comparator<Vertex>{
	@Override
    public int compare(Vertex x, Vertex y)
    {
        if(x.startindex < y.startindex)
        {
        	return -1;
        }
        else if(x.startindex > y.startindex)
        {
        	return 1;
        }
        else
        {
        	if (x.left < y.left)
        	{
        		return -1;
        	}
        	else if (x.left > y.left)
        	{
        		return 1;
        	} 
        	else if(x.left == y.left)
        	{
        		if (x.right < y.right)
            	{
            		return -1;
            	}
            	if (x.right > y.right)
            	{
            		return 1;
            	} 
        		
        	}
        	return 0;
        }
       
    }
}
