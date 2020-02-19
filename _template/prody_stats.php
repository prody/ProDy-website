<?php
$dir = 'sqlite:prody_stats.db';
$dbh  = new PDO($dir) or die("cannot open the database");

$stmt = $dbh->prepare("SELECT SUM(value) AS value_sum FROM codes");
$stmt->execute();

$row = $stmt->fetchAll(PDO::FETCH_OBJ);
$sum = $row->value_sum;

echo $sum;
?>