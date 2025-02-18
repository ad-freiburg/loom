> [!NOTE]
> WORK IN PROGRESS
## Examples

### Filter by route ID

```bash
gtfs2graph --mot 2,bus --route 1,2 examples/zips/sample-feed.zip | topo > sample-feed.json
cat sample-feed.json | loom | transitmap > sample-feed.svg
```