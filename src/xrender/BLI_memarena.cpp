/**
 * $Id: BLI_memarena.c,v 1.5 2005/03/21 00:48:19 zuster Exp $
 *
 * ***** BEGIN GPL/BL DUAL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. The Blender
 * Foundation also sells licenses for use in proprietary software under
 * the Blender License.  See http://www.blender.org/BL/ for information
 * about this.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL/BL DUAL LICENSE BLOCK *****
 * Efficient memory allocation for lots of similar small chunks.
 */
#include "stdafx.h"
#include <stdlib.h>

#include "MEM_guardedalloc.h"

//#include "BLI_blenlib.h"
#include "BLI_memarena.h"
//#include "BLI_linklist.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

struct MemArena {
	unsigned char *curbuf;
	int bufsize, cursize;
	
	LinkNode *bufs;
};

MemArena *BLI_memarena_new(int bufsize) {
	MemArena *ma= (MemArena*)MEM_callocN(sizeof(*ma), "memarena");
	ma->bufsize= bufsize;
	
	return ma;
}
void BLI_memarena_free(MemArena *ma) {
	BLI_linklist_free(ma->bufs, (void(*)(void*)) free);
	MEM_freeN(ma);
}

	/* amt must be power of two */
#define PADUP(num, amt)	((num+(amt-1))&~(amt-1))

void *BLI_memarena_alloc(MemArena *ma, int size) {
	void *ptr;

		/* ensure proper alignment by rounding
		 * size up to multiple of 8 */	
	size= PADUP(size, 8);
	
	if (size>ma->cursize) {
		ma->cursize= (size>ma->bufsize)?size:ma->bufsize;
		ma->curbuf= (unsigned char*)malloc(ma->cursize);
		
		BLI_linklist_prepend(&ma->bufs, ma->curbuf);
	}
	
	ptr= ma->curbuf;
	ma->curbuf+= size;
	ma->cursize-= size;
	
	return ptr;
}

int BLI_linklist_length(LinkNode *list) {
	if (0) {
		return list?(1+BLI_linklist_length(list->next)):0;
	} else {
		int len;

		for (len=0; list; list= list->next)
			len++;

		return len;
	}
}

void BLI_linklist_reverse(LinkNode **listp) {
	LinkNode *rhead= NULL, *cur= *listp;

	while (cur) {
		LinkNode *next= cur->next;

		cur->next= rhead;
		rhead= cur;

		cur= next;
	}

	*listp= rhead;
}

void BLI_linklist_prepend(LinkNode **listp, void *ptr) {
	LinkNode *nlink= (LinkNode*)MEM_mallocN(sizeof(*nlink), "nlink");
	nlink->link= ptr;

	nlink->next= *listp;
	*listp= nlink;
}

void BLI_linklist_append(LinkNode **listp, void *ptr) {
	LinkNode *nlink= (LinkNode*)MEM_mallocN(sizeof(*nlink), "nlink");
	LinkNode *node = *listp;

	nlink->link = ptr;
	nlink->next = NULL;

	if(node == NULL){
		*listp = nlink;
	} else {
		while(node->next != NULL){
			node = node->next;   
		}
		node->next = nlink;
	}
}

void BLI_linklist_prepend_arena(LinkNode **listp, void *ptr, MemArena *ma) {
	LinkNode *nlink= (LinkNode*)BLI_memarena_alloc(ma, sizeof(*nlink));
	nlink->link= ptr;

	nlink->next= *listp;
	*listp= nlink;
}

void BLI_linklist_free(LinkNode *list, LinkNodeFreeFP freefunc) {
	while (list) {
		LinkNode *next= list->next;

		if (freefunc)
			freefunc(list->link);
		MEM_freeN(list);

		list= next;
	}
}

void BLI_linklist_apply(LinkNode *list, LinkNodeApplyFP applyfunc) {
	for (; list; list= list->next)
		applyfunc(list->link);
}



void BLI_addtail(ListBase *listbase, void *vlink)
{
	struct Link *link= (Link*)vlink;

	if (link == 0) return;
	if (listbase == 0) return;

	link->next = 0;
	link->prev = (Link*)listbase->last;

	if (listbase->last) ((struct Link *)listbase->last)->next = link;
	if (listbase->first == 0) listbase->first = link;
	listbase->last = link;
}


void BLI_freelistN(ListBase *listbase)
{
	struct Link *link,*next;

	if (listbase == 0) return;
	link= (Link*)listbase->first;
	while(link) {
		next= link->next;
		MEM_freeN(link);
		link= next;
	}
	listbase->first=0;
	listbase->last=0;
}
