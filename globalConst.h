/*
 * globalConst.h
 *
 * This file is part of the Simon Memory Game project.
 *
 * Copyright 2012 Randal S. Olson.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _globalConst_h_included_
#define _globalConst_h_included_

#define     randDouble      ((double)rand() / (double)RAND_MAX)
#define     maxNodes        256

#define     strLength       6
#define     numInputs       (int)(2 * strLength)
#define     numOutputs      (int)strLength
#define     numColors       2

#define     directedMutations

#endif